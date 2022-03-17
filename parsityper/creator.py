#!/usr/bin/python
import copy
import gzip
import logging
import os
import psutil
import re
import sys
import time
from argparse import (ArgumentParser)
from functools import partial
from mimetypes import guess_type
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from scipy.stats import entropy
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

from parsityper.constants import iupac_replacement, KMER_FIELDS
from parsityper.ext_tools.jellyfish import run_jellyfish_count, parse_jellyfish_counts
from parsityper.helpers import find_overlaping_gene_feature
from parsityper.helpers import init_console_logger, read_tsv, parse_reference_sequence
from parsityper.helpers import read_fasta, generate_non_gap_position_lookup
from parsityper.kmerSearch.kmerSearch import init_automaton_dict, find_in_fasta_dict
from parsityper.version import __version__


def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Parsityper scheme creator')
    parser.add_argument('--input_msa', type=str, required=True,
                        help='MSA in fasta format')
    parser.add_argument('--input_meta', type=str, required=True,
                        help='TSV file of sample_id,genotype')
    parser.add_argument('--ref_id', type=str, required=True,
                        help='sample_id for reference sequence to use in MSA')
    parser.add_argument('--ref_gbk', type=str, required=False,
                        help='GenBank file for reference sequences')
    parser.add_argument('--outdir', type=str, required=True,
                        help='output directory')
    parser.add_argument('--prefix', type=str, required=False,
                        help='output file prefix', default='parsityper')
    parser.add_argument('--min_kmer_freq', type=int, required=False,
                        help='Minimum frequency of kmer for inclusion (default=1)', default=1)
    parser.add_argument('--min_var_freq', type=int, required=False,
                        help='Minimum frequency of variant for inclusion (default=1)', default=1)
    parser.add_argument('--min_ref_frac', type=float, required=False,
                        help='Minimum fraction of isolates positive for reference base for it to be positive 0 - 1.0 (default=1.0)',
                        default=1)
    parser.add_argument('--min_alt_frac', type=float, required=False,
                        help='Minimum fraction of isolates positive for mutation for it to be positive 0 - 1.0 (default=0.95)',
                        default=0.95)
    parser.add_argument('--kmer_len', type=int, required=False,
                        help='Length of kmer to use', default=21)
    parser.add_argument('--max_ambig', type=int, required=False,
                        help='Absolute maximum of degenerate bases allowed in a k-mer (default=0)', default=0)
    parser.add_argument('--max_homo', type=int, required=False,
                        help='Absolute maximum of homopolymer run default=20% of kmer len')
    parser.add_argument('--n_threads', type=int, required=False,
                        help='Number of threads to use', default=1)
    parser.add_argument('--no_plots', required=False,
                        help='suppress making plots, required for large datasets', action='store_true')
    parser.add_argument('--seq_batch_size', type=int, required=False,
                        help='Smaller batch sizes require less memory at the expense of longer run times')
    parser.add_argument('-V', '--version', action='version', version='%(prog)s {}'.format(__version__))
    return parser.parse_args()


def find_gaps(seq):
    """
    Accepts a string and returns the positions of all of the gaps in the sequence
    :param seq: str
    :return: list of [start,end] of all of the gaps
    """
    match = re.finditer(r"\w-+\w", seq)

    positions = []
    for m in match:
        positions.append([m.start() + 1, m.end() - 2])
    return positions

def init_consensus(seq):
    seq_len = len(seq)
    consensus = []
    for i in range(0, seq_len):
        consensus.append({'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0, '-': 0})
    return consensus

def parseVariants(consensus, all_seq_ids, gaps, ref_seq, min_var_freq):
    variants = {'snp': {}, 'del': {}, 'ins': {}}
    bases = ['A', 'T', 'C', 'G']
    for i in range(0, len(ref_seq)):
        pos = consensus[i]
        count_bases = 0
        ref_base = ref_seq[i]
        if ref_base == '-':
            continue
        alt_bases = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
        for b in bases:
            if pos[b] > 0:
                if b != ref_base:
                    alt_bases[b] = pos[b]

        for b in alt_bases:
            if alt_bases[b] < min_var_freq:
                alt_bases[b] = 0
            else:
                count_bases += 1

        if count_bases >= 1:
            variants['snp'][i] = {'ref': ref_base, 'alt': alt_bases}

    for gap in gaps:
        (start, end) = gap.split(':')
        seq_ids = gaps[gap]
        start = int(start)
        end = int(end)
        ref_bases = ref_seq[start:end + 1].replace('-', '')
        vType = 'del'
        if len(ref_bases) == 0:
            vType = 'ins'
            seq_ids = list(set(all_seq_ids) - set(seq_ids))
        if len(seq_ids) < min_var_freq:
            continue
        variants[vType][gap] = {'start': start, 'end': end, 'seq_ids': seq_ids}
    return variants

def preprocess_seqs(fasta_file, ref_id, sample_ids, batch_size, out_dir, prefix, iupac_replacement, min_var_freq):
    trans = str.maketrans(''.join(iupac_replacement.keys()), ''.join(iupac_replacement.values()))
    sample_batch_index = {}
    encoding = guess_type(fasta_file)[1]
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
    is_align_ok = True
    seq_ids = []
    subset_files = []
    file_index = 0
    batch_count = 0
    gaps = {}
    ref_seq = ''
    unalign_files = []
    with _open(fasta_file) as f:
        seq_record = next(SeqIO.parse(f, 'fasta'))
        seq = str(seq_record.seq).upper()
        seq = seq.translate(trans)
        consensus = init_consensus(seq)
        align_len = len(seq)
    out_unalin = []
    out_alin = []
    seq_base_range = range(0, align_len)
    with _open(fasta_file) as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            stime = time.time()
            id = str(seq_record.id)
            if id not in sample_ids:
                continue
            sample_batch_index[id] = file_index
            if batch_count == 0:
                sub_outfile = os.path.join(out_dir, "{}-{}.subset.fasta".format(prefix, file_index))
                unaln_outfile = os.path.join(out_dir, "{}-{}.unaln.fasta".format(prefix, file_index))
                unalign_files.append(unaln_outfile)
                subset_files.append(sub_outfile)
                sFH = open(sub_outfile, 'w')
                uFH = open(unaln_outfile, 'w')
            batch_count += 1
            seq_ids.append(id)
            seq = str(seq_record.seq).upper()
            seq = seq.translate(trans)
            seq_len = len(seq)
            if id == ref_id:
                ref_seq = seq

            intGaps = find_gaps(seq)

            for (start, end) in intGaps:
                gap = "{}:{}".format(start, end)
                if not gap in gaps:
                    gaps[gap] = []
                gaps[gap].append(id)

            for i in seq_base_range:
                consensus[i][seq[i]] += 1

            if seq_len != align_len:
                logging.error("Sequence length mismatch {}: {} vs. {}".format(id, seq_len, align_len))
                is_align_ok = False
            if id not in sample_ids:
                logging.error("Sequence id: {} is not present in metadata file...skip".format(id))
                continue
            out_unalin.append(">{}\n{}".format(id, seq.replace('-', '')))
            out_alin.append(">{}\n{}".format(id, seq))
            if batch_count == batch_size:
                batch_count = 0
                sFH.write("{}\n".format("\n".join(out_alin)))
                uFH.write("{}\n".format("\n".join(out_unalin)))
                sFH.close()
                uFH.close()
                out_unalin = []
                out_alin = []
                file_index += 1
        if len(out_alin) > 0:
            if not sFH.closed:
                sFH.write("{}\n".format("\n".join(out_alin)))
                sFH.close()
        if len(out_unalin) > 0:
            if not uFH.closed:
                uFH.write("{}\n".format("\n".join(out_unalin)))
                uFH.close()

    variants = parseVariants(consensus, seq_ids, gaps, ref_seq, min_var_freq)

    return {'seq_ids': seq_ids, 'is_align_ok': is_align_ok, 'consensus': consensus, 'subset_files': subset_files,
            'unalign_files': unalign_files, 'align_len': align_len,
            'sample_batch_index': sample_batch_index, 'gaps': gaps, 'variants': variants, 'ref_aln': ref_seq}

def perform_kmer_counting(file_manifest, kLen, jellyfish_mem, n_threads):
    if n_threads > 1:
        pool = Pool(processes=n_threads)
    res = []
    if len(file_manifest) < n_threads:
        n_threads = len(file_manifest)
    for i in range(0, len(file_manifest)):
        seq_file = file_manifest[i]
        kmer_file = "{}.jellyfish.txt".format(Path(seq_file).with_suffix(''))
        if n_threads > 1:
            res.append(pool.apply_async(run_jellyfish_count, (seq_file, kmer_file, jellyfish_mem, kLen, n_threads)))
        else:
            res.append(run_jellyfish_count(seq_file, kmer_file, jellyfish_mem, kLen, n_threads))

    if n_threads > 1:
        pool.close()
        pool.join()
        for i in range(0, len(res)):
            res[i].get()
    return res

def combine_jellyfish_results(file_manifest):
    kmers = {}
    for i in range(0, len(file_manifest)):
        kmer_file = file_manifest[i]
        kmer_file = "{}.jellyfish.txt".format(Path(kmer_file).with_suffix(''))
        kmer_df = parse_jellyfish_counts(kmer_file)
        for row in kmer_df.itertuples():
            kmer = row.kmer
            count = row.count
            if not kmer in kmers:
                kmers[kmer] = 0
            kmers[kmer] += count
    return kmers

def calc_homopolymers(seq):
    longest = 0
    for b in ['A', 'T', 'C', 'C']:
        matches = re.findall("{}+".format(b), seq)
        for m in matches:
            length = len(m)
            if length > longest:
                longest = length
    return longest

def filter_kmers(kMers, min_count, max_count, max_ambig, max_homo):
    filtered = {}
    for kmer in kMers:
        count = kMers[kmer]
        if count < min_count or count >= max_count:
            continue
        n_count = kmer.count('N')
        if n_count > max_ambig:
            continue
        homo_len = calc_homopolymers(kmer)
        if homo_len > max_homo:
            continue
        filtered[kmer] = count
    return filtered

def read_fasta(fasta_file):
    """
    Reads fasta file into a dict
    :param fasta_file: fasta formatted file
    :return: dict of sequences
    """
    reference = {}

    encoding = guess_type(fasta_file)[1]
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

    with _open(fasta_file) as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            reference[str(seq_record.id)] = str(seq_record.seq).upper()

    return reference

def write_aho_kmerResults(aho, seq_file, out_file):
    find_in_fasta_dict(aho, read_fasta(seq_file)).reset_index().to_csv(out_file, header=True, sep="\t")

def process_aho_kmerResults(aho, num_kmers, seq_file, kLen):
    seqs = read_fasta(seq_file)
    kmer_results = {}
    sample_results = {}
    iter_keys = seqs.keys()
    aln_lookup = {}
    for seq_id in iter_keys:
        sample_results[seq_id] = [0] * num_kmers
        seq = seqs[seq_id].replace('-', '')
        for idx, (kIndex, kmer_seq, is_revcomp) in aho.iter(seq):
            kIndex = int(kIndex)
            sample_results[seq_id][kIndex] = 1
            if not kIndex in kmer_results:
                kmer_results[kIndex] = {'kSeq': kmer_seq, 'aSeq': '', 'rSeq': '', 'match_index': idx,
                                        'is_revcomp': is_revcomp, 'uStart': 0, 'uEnd': 0, 'aStart': 0, 'aEnd': 0,
                                        'seq_id': seq_id}
                if not seq_id in aln_lookup:
                    aln_lookup[seq_id] = create_aln_pos_from_unalign_pos(seqs[seq_id])

    for kIndex in kmer_results:
        seq_id = kmer_results[kIndex]['seq_id']
        match_index = int(kmer_results[kIndex]['match_index'])
        unalign_pos_end = match_index
        unalign_pos_start = match_index - kLen + 1

        align_pos_start = aln_lookup[seq_id][unalign_pos_start]
        align_pos_end = aln_lookup[seq_id][unalign_pos_end]
        aSeq = seqs[seq_id][align_pos_start:align_pos_end + 1]
        kmer_results[kIndex]['aSeq'] = aSeq
        kmer_results[kIndex]['uStart'] = unalign_pos_start
        kmer_results[kIndex]['uEnd'] = unalign_pos_end
        kmer_results[kIndex]['aStart'] = align_pos_start
        kmer_results[kIndex]['aEnd'] = align_pos_end

    return {'kmer': kmer_results, 'samples': sample_results}

def map_kmers(kMers, file_manifest, n_threads):
    kmer_index = {}
    i = 0
    for kmer in kMers:
        kmer_index[i] = kmer
        i += 1
    kLen = len(kmer_index[0])
    num_kmers = len(kmer_index)
    aho = init_automaton_dict(kmer_index)

    if n_threads > 1:
        pool = Pool(processes=n_threads)
    results = []
    for i in range(0, len(file_manifest)):
        seq_file = file_manifest[i]
        if n_threads > 1:
            results.append(pool.apply_async(process_aho_kmerResults, (aho, num_kmers, seq_file, kLen)))
        else:
            results.append(process_aho_kmerResults(aho, num_kmers, seq_file, kLen))

    if n_threads > 1:
        pool.close()
        pool.join()

        for i in range(0, len(results)):
            results[i] = results[i].get()

    return combine_aho_results(results)

def combine_aho_results(results):
    kmers = {}
    samples = {}
    for i in range(0, len(results)):
        for kmer in results[i]['kmer']:
            if kmer in kmers:
                continue
            kmers[kmer] = results[i]['kmer'][kmer]
        for sample_id in results[i]['samples']:
            samples[sample_id] = results[i]['samples'][sample_id]
    return {'kmer': kmers, 'samples': samples}

def add_ref_kmer_info(kmer_results, kmer_counts, ref_seq, kLen):
    kmer_base_range = range(0, kLen)
    for i in kmer_results['kmer']:
        aStart = kmer_results['kmer'][i]['aStart']
        aEnd = kmer_results['kmer'][i]['aEnd']
        kmer_results['kmer'][i]['rSeq'] = ref_seq[aStart:aEnd + 1]
        count_var = 0
        for k in kmer_base_range:
            if kmer_results['kmer'][i]['rSeq'][k] != kmer_results['kmer'][i]['aSeq'][k]:
                count_var += 1
        kmer_results['kmer'][i]['count_var'] = count_var
        kmer_results['kmer'][i]['total_count'] = kmer_counts[kmer_results['kmer'][i]['kSeq']]
        kmer_results['kmer'][i]['positive_genotypes'] = []
        kmer_results['kmer'][i]['partial_genotypes'] = []

def get_kmer_genotype_counts_bck(kmer_files, genotypeMapping):
    genotypes = list(set(genotypeMapping.values()))
    num_genotypes = len(genotypes)
    genotype_kCounts = {}
    for i in range(0, len(kmer_files)):
        kmer_df = read_tsv(kmer_files[i])
        for row in kmer_df.itertuples():
            kmer = row.kmername
            seq_id = row.contig_id
            pos = row.match_index
            genotype = genotypeMapping[seq_id]
            if not kmer in genotype_kCounts:
                genotype_kCounts[kmer] = {'total': 0, 'counts': {}, 'pos': pos, 'seq_id': seq_id}
                for k in range(0, num_genotypes):
                    genotype_kCounts[kmer]['counts'][genotypes[k]] = 0
            genotype_kCounts[kmer]['counts'][genotype] += 1
            genotype_kCounts[kmer]['total'] += 1
    return genotype_kCounts

def get_kmer_genotype_counts(valid_kmer_indicies, sample_profiles, genotypeMapping):
    genotypes = list(set(genotypeMapping.values()))
    genotype_kCounts = {}
    for i in valid_kmer_indicies:
        genotype_kCounts[i] = {'total': 0, 'counts': {}}
        for genotype in genotypes:
            genotype_kCounts[i]['counts'][genotype] = 0

    num_kmers = len(valid_kmer_indicies)
    kmer_range = range(0, num_kmers)
    for sample_id in sample_profiles:
        genotype = genotypeMapping[sample_id]
        for i in kmer_range:
            index = valid_kmer_indicies[i]
            value = sample_profiles[sample_id][index]
            genotype_kCounts[index]['counts'][genotype] += value
            genotype_kCounts[index]['total'] += value
    return genotype_kCounts

def calc_shanon_entropy(value_list):
    total = sum(value_list)
    if total == 0:
        return -1
    values = []
    for v in value_list:
        values.append(v / total)
    return entropy(values)

def calc_AMI(category_1, category_2):
    return adjusted_mutual_info_score(category_1, category_2, average_method='arithmetic')

def calc_ARI(category_1, category_2):
    return adjusted_rand_score(category_1, category_2)

def calc_kmer_entropy(genotype_kCounts):
    sEntropy = {}
    for kmer in genotype_kCounts:
        counts = list(genotype_kCounts[kmer]['counts'].values())
        sEntropy[kmer] = calc_shanon_entropy(counts)
    return sEntropy

def calc_kmer_associations_bck(genotype_kCounts, genotypeCounts, sample_count, n_threads=1):
    sampleInfo = {}
    sample_padding = sample_count * (len(genotypeCounts) - 1)

    for kmer in genotype_kCounts:
        sampleInfo[kmer] = {}
        kVec = []
        genotypeVec = {}

        # Init the kmer presence vector accross all samples
        for genotype in genotype_kCounts[kmer]['counts']:
            total = genotypeCounts[genotype]
            num_pos = total - genotype_kCounts[kmer]['counts'][genotype]
            perc_pos = num_pos / total
            scaled_num_pos = int(perc_pos * sample_count)
            scaled_num_neg = sample_count - scaled_num_pos
            gVec = ([1] * scaled_num_pos) + ([0] * scaled_num_neg)
            kVec += gVec
            genotypeVec[genotype] = gVec + [0] * sample_padding

        # Compare ARI and AMI for each kmer with each genotype
        for genotype in genotypeVec:
            gVec = genotypeVec[genotype]
            if len(gVec) != len(kVec):
                print("{}\n{}\n{}\n".format(kmer, kVec, gVec))
            ami = calc_AMI(gVec, kVec)
            ari = calc_ARI(gVec, kVec)
            sampleInfo[kmer][genotype] = {'ari': ari, 'ami': ami}
    return sampleInfo

def calc_kmer_associations(query_genotype, kmer, genotype_kCounts, genotypeCounts, sample_count):
    sample_padding = sample_count * (len(genotypeCounts) - 1)
    kVec = []
    genotypeVec = {}
    # Init the kmer presence vector accross all samples
    for genotype in genotype_kCounts[kmer]['counts']:
        total = genotypeCounts[genotype]
        num_pos = total - genotype_kCounts[kmer]['counts'][genotype]
        perc_pos = num_pos / total
        scaled_num_pos = int(perc_pos * sample_count)
        scaled_num_neg = sample_count - scaled_num_pos
        gVec = ([1] * scaled_num_pos) + ([0] * scaled_num_neg)
        kVec += gVec
        genotypeVec[genotype] = gVec + [0] * sample_padding

    # Compare ARI and AMI for each kmer with each genotype

    gVec = genotypeVec[query_genotype]
    if len(gVec) != len(kVec):
        logging.error("{}\n{}\n{}\n".format(kmer, kVec, gVec))
    ami = calc_AMI(gVec, kVec)
    ari = calc_ARI(gVec, kVec)

    return {'ami':ami,'ari':ari}

def get_genotype_mapping(metadata_df):
    '''

    :param metadata_df:
    :type metadata_df:
    :return:
    :rtype:
    '''
    mapping = {}
    for row in metadata_df.itertuples():
        sample_id = row.sample_id
        genotype = row.genotype
        mapping[sample_id] = genotype
    return mapping

def get_non_gap_position(ref_non_gap_lookup, pos):
    non_gap_position = ref_non_gap_lookup[pos]
    while non_gap_position == -1:
        pos -= 1
        non_gap_position = ref_non_gap_lookup[pos]
    return non_gap_position

def create_aln_pos_from_unalign_pos(aln_seq):
    unalign_seq = aln_seq.replace('-', '')
    aln_len = len(aln_seq)
    unaln_len = len(unalign_seq)
    lookup = [-1] * unaln_len
    pos = 0
    for i in range(0, unaln_len):
        for k in range(pos, aln_len):
            if unalign_seq[i] == aln_seq[k]:
                lookup[i] = k
                pos = k + 1
                break
    return lookup

def generate_non_gap_position_lookup(seq):
    """
    Creates a list of positions which correspond to the position of that base in a gapless sequence
    :param seq: string
    :return: list
    """
    length = len(seq)
    num_gaps = 0
    lookup = []
    for i in range(0, length):
        base = seq[i]
        if base == '-':
            num_gaps += 1
            lookup.append(-1)
        else:
            lookup.append(i - num_gaps)
    return lookup

def get_kmer_positions(ref_id, genotype_kCounts, kLen, fasta_dir, individual_file_sample_index):
    positions = {}
    aln_lookup = {}
    non_gap_lookup = {}
    file_index = individual_file_sample_index[ref_id]
    fasta_file = os.path.join(fasta_dir, "{}/{}.fasta".format(file_index, ref_id))
    ref_seq = read_fasta(fasta_file)[ref_id]
    aln_len = len(ref_seq)
    aln_lookup[ref_id] = create_aln_pos_from_unalign_pos(ref_seq)
    non_gap_lookup[ref_id] = generate_non_gap_position_lookup(ref_seq)

    # saves having to sort after the fact
    for i in range(0, aln_len):
        positions[i] = {}
    seqs = {ref_id: ref_seq}

    for kmer in genotype_kCounts:
        unalign_pos_start = genotype_kCounts[kmer]['pos'] - kLen + 1
        unalign_pos_end = unalign_pos_start + kLen - 1
        seq_id = genotype_kCounts[kmer]['seq_id']
        file_index = individual_file_sample_index[seq_id]
        if not seq_id in aln_lookup:
            fasta_file = os.path.join(fasta_dir, "{}/{}.fasta".format(file_index, seq_id))
            seqs[seq_id] = read_fasta(fasta_file)[seq_id]
            aln_lookup[seq_id] = create_aln_pos_from_unalign_pos(seqs[seq_id])
        align_pos_start = aln_lookup[seq_id][unalign_pos_start]
        align_pos_end = aln_lookup[seq_id][unalign_pos_end]
        kTest = seqs[seq_id].replace('-', '')[unalign_pos_start:unalign_pos_end + 1]
        is_rev = False
        if kTest != kmer:
            is_rev = True
        align_refSeq_kmer = ref_seq[align_pos_start:align_pos_end + 1]
        align_querySeq_kmer = seqs[seq_id][align_pos_start:align_pos_end + 1]
        kmer_bases = []
        count_var = 0

        for i in range(0, kLen):
            if align_refSeq_kmer[i] != align_querySeq_kmer[i]:
                if align_refSeq_kmer[i] != 'N' and align_querySeq_kmer[i] != 'N':
                    count_var += 1

        unalign_pos_start = get_non_gap_position(non_gap_lookup[ref_id], unalign_pos_start)
        unalign_pos_end = get_non_gap_position(non_gap_lookup[ref_id], unalign_pos_end)
        positions[align_pos_start][kmer] = {'aStart': align_pos_start, 'aEnd': align_pos_end,
                                            'uStart': unalign_pos_start,
                                            'uEnd': unalign_pos_end, 'count_var_pos': count_var, 'is_rev': is_rev,
                                            'entropy': -1, 'total_count': 0,
                                            'align_refSeq_kmer': align_refSeq_kmer,
                                            'align_querySeq_kmer': align_querySeq_kmer, 'seq_id': seq_id,
                                            'positive_genotypes': [], 'partial_genotypes': []}
    return positions

def add_gene_inference(selected_kmers, ref_seq, ref_id, reference_info, trans_table=1):

    aln_refLookup = create_aln_pos_from_unalign_pos(ref_seq)
    for mutation_type in selected_kmers:
        for event in selected_kmers[mutation_type]:
            ref_var = ''
            alt_var = ''
            aa_name = ''
            gene_name = 'intergenic'
            gene_len = 0
            gene_start = -1
            gene_end = -1
            gene_seq = ''
            positions = ''
            is_cds = False
            is_silent = True
            is_frameshift = False

            if mutation_type == 'snp':
                start = event
                end = event
                for base in selected_kmers[mutation_type][event]['ref']['kmers']:
                    if len(selected_kmers[mutation_type][event]['ref']['kmers'][base]) > 0:
                        ref_var = base
                var_len = 1
            else:
                (start, end) = event.split(':')
                start = int(start)
                end = int(end)
                ref_var = ref_seq[start:end + 1]
                var_len = len(ref_var.replace('-', ''))

            ref_var_dna = ref_var
            gene_feature = find_overlaping_gene_feature(start, end, reference_info, ref_id)
            ref_var_aa = ''
            if gene_feature is not None:
                is_cds = True
                if mutation_type != 'snp' and (var_len - 1) % 3 > 0:
                    is_frameshift = True

                gene_name = gene_feature['gene_name']
                gene_len = gene_feature['gene_len']
                positions = gene_feature['positions']
                for s, e in positions:
                    if gene_start == -1:
                        gene_start = s
                        gene_end = e
                    else:
                        if gene_start > s:
                            gene_start = s
                        if gene_end < e:
                            gene_end = e
                aln_gene_start = aln_refLookup[gene_start]
                aln_gene_end = aln_refLookup[gene_end]
                if aln_gene_start == aln_gene_end:
                    aln_gene_seq = ref_seq[aln_gene_start:aln_gene_end + 1]
                else:
                    aln_gene_seq = ref_seq[aln_gene_start:aln_gene_end]
                rel_var_start = abs(start - aln_gene_start)
                rel_var_end = abs(end - aln_gene_start)
                aa_var_start = int(rel_var_start / 3)
                codon_var_start = aa_var_start * 3
                var_len = rel_var_end - rel_var_start + 1
                rem = 3 - var_len % 3
                codon_var_end = int(codon_var_start + (var_len + rem))
                aa_var_end = int(codon_var_end / 3)
                if codon_var_end == 0:
                    codon_var_end = codon_var_start + 3
                ref_var_dna = aln_gene_seq[codon_var_start:codon_var_end].replace('-', '')
                incr=0
                while len(ref_var_dna) % 3 != 0:
                    incr += 1
                    ref_var_dna = ''.join(aln_gene_seq[codon_var_start:codon_var_end + incr]).replace('-', '')
                    if incr + codon_var_end > gene_len:
                        break
                ref_var_aa = str(Seq(ref_var_dna).translate(table=trans_table))

            if mutation_type == 'snp':
                for kmer in selected_kmers[mutation_type][event]['ref']['kmers'][ref_var]:

                    selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['is_cds'] = is_cds
                    selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['is_silent'] = is_silent
                    selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['is_frameshift'] = is_frameshift
                    selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['gene_name'] = gene_name
                    selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['gene_len'] = gene_len
                    selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['ref_var'] = ref_var
                    selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['target_var'] = ref_var
                    selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['ref_var_dna'] = ref_var_dna
                    selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['ref_var_aa'] = ref_var_aa
                    if is_cds:
                        selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer][
                            'cds_start'] = codon_var_start
                        selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer][
                            'cds_end'] = codon_var_end - 1
                        selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['aa_start'] = aa_var_start
                        selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['aa_end'] = aa_var_end - 1
                        selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer][
                            'aa_name'] = "{}{}{}".format(ref_var_aa, aa_var_start, ref_var_aa)
            else:
                for kmer in selected_kmers[mutation_type][event]['ref']['kmers']:
                    selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['is_cds'] = is_cds
                    selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['is_silent'] = is_silent
                    selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['is_frameshift'] = is_frameshift
                    selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['gene_name'] = gene_name
                    selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['gene_len'] = gene_len
                    selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['ref_var'] = ref_var
                    selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['target_var'] = ref_var
                    selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['ref_var_dna'] = ref_var_dna
                    selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['ref_var_aa'] = ref_var_aa
                    if is_cds:
                        selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['cds_start'] = codon_var_start
                        selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['cds_end'] = codon_var_end - 1
                        selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['aa_start'] = aa_var_start
                        selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['aa_end'] = aa_var_end - 1
                        selected_kmers[mutation_type][event]['ref']['kmers'][kmer][
                            'aa_name'] = "{}_{}_{}".format(mutation_type, aa_var_start, aa_var_end - 1)

            if mutation_type == 'snp':
                for base in selected_kmers[mutation_type][event]['alt']['kmers']:
                    if len(selected_kmers[mutation_type][event]['alt']['kmers'][base]) == 0:
                        continue
                    for kmer in selected_kmers[mutation_type][event]['alt']['kmers'][base]:
                        align_querySeq_kmer = selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['aSeq']
                        akmer_start = selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['aStart']
                        akmer_end = selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['aEnd']
                        rel_alt_start = abs(start - akmer_start)
                        rel_alt_end = rel_alt_start + var_len
                        alt_var = align_querySeq_kmer[rel_alt_start:rel_alt_end]
                        rel_alt_start = abs(start - gene_start)
                        alt_var_aa = ''
                        alt_var_dna = alt_var
                        alt_seq = list(ref_seq)


                        if is_cds:
                            for i in range(akmer_start, akmer_end + 1):
                                pos = i - akmer_start
                                alt_seq[i] = align_querySeq_kmer[pos]

                            alt_gene_seq = ''.join(alt_seq)[gene_start:]
                            alt_gene_seq = alt_gene_seq.replace('-','')[0:gene_len]
                            alt_gene_aa = str(Seq(alt_gene_seq).translate(table=trans_table))
                            alt_var_dna = ''.join(alt_gene_seq[codon_var_start:codon_var_end]).replace('-', '')
                            alt_var_aa = alt_gene_aa[aa_var_start:aa_var_end]
                            selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['cds_start'] = codon_var_start
                            selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['cds_end'] = codon_var_end - 1
                            selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['aa_start'] = aa_var_start
                            selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['aa_end'] = aa_var_end - 1
                            selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['aa_name'] = "{}{}{}".format(ref_var_aa, aa_var_start, alt_var_aa)
                            if alt_var_aa != ref_var_aa:
                                is_silent = False
                        selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['is_silent'] = is_silent
                        selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['is_frameshift'] = is_frameshift
                        selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['is_cds'] = is_cds
                        selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['gene_name'] = gene_name
                        selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['gene_len'] = gene_len
                        selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['ref_var'] = ref_var
                        selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['target_var'] = alt_var
                        selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['alt_var_dna'] = alt_var_dna
                        selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['alt_var_aa'] = alt_var_aa
            else:
                for kmer in selected_kmers[mutation_type][event]['alt']['kmers']:
                    align_querySeq_kmer = selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['aSeq']
                    akmer_start = selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['aStart']
                    akmer_end = selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['aEnd']
                    rel_alt_start = abs(start - akmer_start)
                    rel_alt_end = rel_alt_start + var_len
                    alt_var = align_querySeq_kmer[rel_alt_start:rel_alt_end]
                    alt_var_aa = ''
                    alt_var_dna = alt_var
                    if is_cds:
                        alt_seq = list(ref_seq)
                        for i in range(akmer_start, akmer_end + 1):
                            pos = i - akmer_start
                            alt_seq[i] = align_querySeq_kmer[pos]

                        alt_gene_seq = ''.join(alt_seq)[gene_start:]
                        alt_gene_seq = alt_gene_seq.replace('-', '')[0:gene_len]
                        alt_gene_aa = str(Seq(alt_gene_seq).translate(table=trans_table))
                        alt_var_dna = ''.join(alt_gene_seq[codon_var_start:codon_var_end]).replace('-', '')
                        alt_var_aa = alt_gene_aa[aa_var_start:aa_var_end]
                        selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['cds_start'] = codon_var_start
                        selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['cds_end'] = codon_var_end - 1
                        selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['aa_start'] = aa_var_start
                        selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['aa_end'] = aa_var_end - 1
                        selected_kmers[mutation_type][event]['alt']['kmers'][kmer][
                            'aa_name'] = "{}_{}_{}".format(mutation_type, aa_var_start, aa_var_end - 1)
                        if alt_var_aa != ref_var_aa:
                            is_silent = False

                    selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['is_silent'] = is_silent
                    selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['is_frameshift'] = is_frameshift
                    selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['is_cds'] = is_cds
                    selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['gene_name'] = gene_name
                    selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['gene_len'] = gene_len
                    selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['ref_var'] = ref_var
                    selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['target_var'] = alt_var
                    selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['alt_var_dna'] = alt_var_dna
                    selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['alt_var_aa'] = alt_var_aa

    return selected_kmers

def add_gene_inference_bck(selected_kmers, ref_seq, ref_id, reference_info, trans_table=1):

    aln_refLookup = create_aln_pos_from_unalign_pos(ref_seq)
    for mutation_type in selected_kmers:
        for event in selected_kmers[mutation_type]:
            ref_var = ''
            alt_var = ''
            aa_name = ''
            gene_name = 'intergenic'
            gene_len = 0
            gene_start = -1
            gene_end = -1
            gene_seq = ''
            positions = ''
            is_cds = False
            is_silent = True
            is_frameshift = False

            if mutation_type == 'snp':
                start = event
                end = event
                for base in selected_kmers[mutation_type][event]['ref']['kmers']:
                    if len(selected_kmers[mutation_type][event]['ref']['kmers'][base]) > 0:
                        ref_var = base
                var_len = 1
            else:
                (start, end) = event.split(':')
                start = int(start)
                end = int(end)
                ref_var = ref_seq[start:end + 1]
                var_len = len(ref_var.replace('-', ''))

            ref_var_dna = ref_var
            gene_feature = find_overlaping_gene_feature(start, end, reference_info, ref_id)
            ref_var_aa = ''
            if gene_feature is not None:
                is_cds = True
                if mutation_type != 'snp' and (var_len - 1) % 3 > 0:
                    is_frameshift = True

                gene_name = gene_feature['gene_name']
                gene_len = gene_feature['gene_len']
                positions = gene_feature['positions']
                for s, e in positions:
                    if gene_start == -1:
                        gene_start = s
                        gene_end = e
                    else:
                        if gene_start > s:
                            gene_start = s
                        if gene_end < e:
                            gene_end = e
                aln_gene_start = aln_refLookup[gene_start]
                aln_gene_end = aln_refLookup[gene_end]
                if aln_gene_start == aln_gene_end:
                    aln_gene_seq = ref_seq[aln_gene_start:aln_gene_end + 1]
                else:
                    aln_gene_seq = ref_seq[aln_gene_start:aln_gene_end]
                rel_var_start = abs(start - aln_gene_start)
                rel_var_end = abs(end - aln_gene_start)
                aa_var_start = int(rel_var_start / 3)
                codon_var_start = aa_var_start * 3
                var_len = rel_var_end - rel_var_start + 1
                rem = 3 - var_len % 3
                codon_var_end = int(codon_var_start + (var_len + rem))
                aa_var_end = int(codon_var_end / 3)
                if codon_var_end == 0:
                    codon_var_end = codon_var_start + 3
                ref_var_dna = aln_gene_seq[codon_var_start:codon_var_end].replace('-', '')
                incr=0
                while len(ref_var_dna) % 3 != 0:
                    incr += 1
                    ref_var_dna = ''.join(aln_gene_seq[codon_var_start:codon_var_end + incr]).replace('-', '')
                    if incr + codon_var_end > gene_len:
                        break
                ref_var_aa = str(Seq(ref_var_dna).translate(table=trans_table))

            if mutation_type == 'snp':
                for kmer in selected_kmers[mutation_type][event]['ref']['kmers'][ref_var]:

                    selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['is_cds'] = is_cds
                    selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['is_silent'] = is_silent
                    selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['is_frameshift'] = is_frameshift
                    selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['gene_name'] = gene_name
                    selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['gene_len'] = gene_len
                    selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['ref_var'] = ref_var
                    selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['target_var'] = ref_var
                    selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['ref_var_dna'] = ref_var_dna
                    selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['ref_var_aa'] = ref_var_aa
                    if is_cds:
                        selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer][
                            'cds_start'] = codon_var_start
                        selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer][
                            'cds_end'] = codon_var_end - 1
                        selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['aa_start'] = aa_var_start
                        selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer]['aa_end'] = aa_var_end - 1
                        selected_kmers[mutation_type][event]['ref']['kmers'][ref_var][kmer][
                            'aa_name'] = "{}{}{}".format(ref_var_aa, aa_var_start, ref_var_aa)
            else:
                for kmer in selected_kmers[mutation_type][event]['ref']['kmers']:
                    selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['is_cds'] = is_cds
                    selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['is_silent'] = is_silent
                    selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['is_frameshift'] = is_frameshift
                    selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['gene_name'] = gene_name
                    selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['gene_len'] = gene_len
                    selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['ref_var'] = ref_var
                    selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['target_var'] = ref_var
                    selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['ref_var_dna'] = ref_var_dna
                    selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['ref_var_aa'] = ref_var_aa
                    if is_cds:
                        selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['cds_start'] = codon_var_start
                        selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['cds_end'] = codon_var_end - 1
                        selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['aa_start'] = aa_var_start
                        selected_kmers[mutation_type][event]['ref']['kmers'][kmer]['aa_end'] = aa_var_end - 1
                        selected_kmers[mutation_type][event]['ref']['kmers'][kmer][
                            'aa_name'] = "{}_{}_{}".format(mutation_type, aa_var_start, aa_var_end - 1)

            if mutation_type == 'snp':
                for base in selected_kmers[mutation_type][event]['alt']['kmers']:
                    if len(selected_kmers[mutation_type][event]['alt']['kmers'][base]) == 0:
                        continue
                    for kmer in selected_kmers[mutation_type][event]['alt']['kmers'][base]:
                        align_querySeq_kmer = selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['aSeq']
                        akmer_start = selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['aStart']
                        akmer_end = selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['aEnd']
                        rel_alt_start = abs(start - akmer_start)
                        rel_alt_end = rel_alt_start + var_len
                        alt_var = align_querySeq_kmer[rel_alt_start:rel_alt_end]
                        rel_alt_start = abs(start - gene_start)
                        alt_var_aa = ''
                        alt_var_dna = alt_var
                        if is_cds:
                            alt_seq = list(aln_gene_seq)

                            for i in range(0, len(alt_var_dna)):
                                pos = i + rel_alt_start

                                alt_seq[pos] = alt_var_dna[i]
                                print("{}\t{}\t{}\t{}\t{}".format(pos, i, alt_var_dna, len(alt_seq), ''.join(alt_seq)))
                            alt_var_dna = ''.join(alt_seq[codon_var_start:codon_var_end]).replace('-', '')
                            alt_var_aa = str(Seq(alt_var_dna).translate(table=trans_table))
                            selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer][
                                'cds_start'] = codon_var_start
                            selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer][
                                'cds_end'] = codon_var_end - 1
                            selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['aa_start'] = aa_var_start
                            selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['aa_end'] = aa_var_end - 1
                            selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer][
                                'aa_name'] = "{}{}{}".format(ref_var_aa, aa_var_start, alt_var_aa)
                            if alt_var_aa != ref_var_aa:
                                is_silent = False
                        selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['is_silent'] = is_silent
                        selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer][
                            'is_frameshift'] = is_frameshift
                        selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['is_cds'] = is_cds
                        selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['gene_name'] = gene_name
                        selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['gene_len'] = gene_len
                        selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['ref_var'] = ref_var
                        selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['target_var'] = alt_var
                        selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['alt_var_dna'] = alt_var_dna
                        selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['alt_var_aa'] = alt_var_aa
            else:
                for kmer in selected_kmers[mutation_type][event]['alt']['kmers']:
                    align_querySeq_kmer = selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['aSeq']
                    akmer_start = selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['aStart']
                    akmer_end = selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['aEnd']
                    rel_alt_start = abs(start - akmer_start)
                    rel_alt_end = rel_alt_start + var_len
                    alt_var = align_querySeq_kmer[rel_alt_start:rel_alt_end]
                    alt_var_aa = ''
                    alt_var_dna = alt_var
                    if is_cds:
                        alt_seq = list(aln_gene_seq)
                        #print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(gene_name,alt_seq,len(alt_seq),start,end,rel_alt_start,rel_alt_end))
                        rel_alt_start = abs(start - gene_start)
                        for i in range(0, len(alt_var_dna)):
                            pos = i + rel_alt_start
                            alt_seq[pos] = alt_var_dna[i]
                        alt_var_dna = ''.join(alt_seq[codon_var_start:codon_var_end]).replace('-', '')
                        incr = 0
                        while len(alt_var_dna) % 3 != 0:
                            incr += 1
                            alt_var_dna = ''.join(alt_seq[codon_var_start:codon_var_end + incr]).replace('-', '')
                            if incr + codon_var_end > gene_len:
                                #print("{}\t{}\t{}\t{}".format(start,end,ref_var_dna,alt_var_dna))
                                break
                        alt_var_aa = str(Seq(alt_var_dna.replace('-', '')).translate(table=trans_table))
                        selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['cds_start'] = codon_var_start
                        selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['cds_end'] = codon_var_end - 1
                        selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['aa_start'] = aa_var_start
                        selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['aa_end'] = aa_var_end - 1
                        selected_kmers[mutation_type][event]['alt']['kmers'][kmer][
                            'aa_name'] = "{}_{}_{}".format(mutation_type, aa_var_start, aa_var_end - 1)
                        if alt_var_aa != ref_var_aa:
                            is_silent = False

                    selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['is_silent'] = is_silent
                    selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['is_frameshift'] = is_frameshift
                    selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['is_cds'] = is_cds
                    selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['gene_name'] = gene_name
                    selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['gene_len'] = gene_len
                    selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['ref_var'] = ref_var
                    selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['target_var'] = alt_var
                    selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['alt_var_dna'] = alt_var_dna
                    selected_kmers[mutation_type][event]['alt']['kmers'][kmer]['alt_var_aa'] = alt_var_aa

    return selected_kmers


def select_snp_kmers(snp_variants, kLen, kmer_alignment_positions):
    variant_kmers = {}
    for aln_start in kmer_alignment_positions:
        kmers = kmer_alignment_positions[aln_start]
        affected_bases = range(aln_start, aln_start + kLen)
        ovlVariants = {}
        for i in affected_bases:
            if i in snp_variants:
                ovlVariants[i] = snp_variants[i]
        if len(ovlVariants) == 0:
            continue

        max_count_var = 0

        for kmer in kmers:
            count_variants = kmers[kmer]['count_var']
            if count_variants > max_count_var:
                max_count_var = count_variants

        for pos in ovlVariants:
            is_new = False
            if not pos in variant_kmers:
                variant_kmers[pos] = {
                    'ref': {'count_var': 0, 'kmers': {'A': {}, 'T': {}, 'C': {}, 'G': {}}},
                    'alt': {'count_var': 0, 'kmers': {'A': {}, 'T': {}, 'C': {}, 'G': {}}}
                }
                #is_new = True
            else:
                #prev_count = variant_kmers[pos]['alt']['count_var']
                #if prev_count > max_count_var:
                continue
                #else:
                #    variant_kmers[pos][state]['kmers'][b] = {}
            ref_base = snp_variants[pos]['ref']
            if ref_base == '-':
                continue
            alt_bases = snp_variants[pos]['alt']
            rel_pos = pos - aln_start

            for kmer in kmers:
                seq = kmers[kmer]['aSeq']
                b = seq[rel_pos]
                if b == '-':
                    continue
                state = 'ref'
                if b != 'N' and b != ref_base and alt_bases[b] > 0:
                    state = 'alt'

                variant_kmers[pos][state]['count_var'] = max_count_var
                variant_kmers[pos][state]['kmers'][b][int(kmer)] = copy.deepcopy(kmer_alignment_positions[aln_start][kmer])


    return variant_kmers

def select_indel_kmers(indel_variants, kLen, kmer_alignment_positions):
    variant_kmers = {}
    max_key = max(list(kmer_alignment_positions.keys()))
    for indel in indel_variants:
        aStart = indel_variants[indel]['start']
        s = aStart - kLen - 1
        if s < 0:
            s = 0
        aEnd = indel_variants[indel]['end']
        e = aEnd + kLen
        if e > max_key:
            e = max_key + 1

        affected_bases = range(s, e)
        variant_kmers[indel] = {'ref': {'count_var': 0, 'kmers': {}}, 'alt': {'count_var': 0, 'kmers': {}}}
        indel_len = aEnd - aStart + 1

        for i in affected_bases:
            if not i in kmer_alignment_positions:
                continue

            if len(kmer_alignment_positions[i]) == 1:
                continue

            kmers = kmer_alignment_positions[i]
            rel_start = abs(aStart - i)
            rel_end = rel_start + indel_len
            max_count_var = 0
            to_add = False
            rKmers = []
            aKmers = []
            for kmer in kmers:
                count_variants = kmers[kmer]['count_var']
                rSeq = kmers[kmer]['rSeq'][rel_start:rel_end].replace('-', '')
                qSeq = kmers[kmer]['aSeq'][rel_start:rel_end].replace('-', '')
                qLen = len(qSeq)
                rLen = len(rSeq)
                delta = abs(rLen - qLen)

                if delta == 0:
                    rKmers.append(kmer)
                elif delta == indel_len:
                    aKmers.append(kmer)
                else:
                    rKmers.append(kmer)

                if count_variants > max_count_var:
                    max_count_var = count_variants
                    to_add = True

            if len(rKmers) == 0 or len(aKmers) == 0:
                continue

            if to_add:
                variant_kmers[indel] = {'ref': {'count_var': 0, 'kmers': {}}, 'alt': {'count_var': 0, 'kmers': {}}}
                variant_kmers[indel]['ref']['count_var'] = max_count_var
                variant_kmers[indel]['alt']['count_var'] = max_count_var
                for kmer in rKmers:
                    variant_kmers[indel]['ref']['kmers'][int(kmer)] = copy.deepcopy(kmer_alignment_positions[i][kmer])
                for kmer in aKmers:
                    variant_kmers[indel]['alt']['kmers'][int(kmer)] = copy.deepcopy(kmer_alignment_positions[i][kmer])

    return variant_kmers

def getSelectedKmerIndicies(selected_kmers):
    indicies = []
    for vType in selected_kmers:
        for event in selected_kmers[vType]:
            for state in ['ref', 'alt']:
                if vType == 'snp':
                    for base in selected_kmers[vType][event][state]["kmers"]:
                        indicies += list(selected_kmers[vType][event][state]["kmers"][base].keys())
                else:
                    indicies += list(selected_kmers[vType][event][state]["kmers"])
    return sorted(list(set(indicies)))

def populate_fields(uid, mutation_type, vStart, vEnd, state, ref_variant, alt_variant, ref_aa,target_aa,kmer_info):
    record = copy.deepcopy(KMER_FIELDS)
    vStart+=1
    vEnd+=1
    for field in record:
        if field in kmer_info:
            record[field] = kmer_info[field]
    record['key'] = uid
    record['mutation_type'] = mutation_type
    record['state'] = state
    record['ref_state'] = ref_variant
    record['alt_state'] = alt_variant
    record['variant_start'] = vStart
    record['variant_end'] = vEnd
    record['kmer_start'] = kmer_info['aStart']
    record['kseq'] = kmer_info['aSeq'].replace('-', '')
    record['klen'] = len(record['kseq'])
    record['homopolymer_len'] = calc_homopolymers(record['kseq'])
    record['target_variant'] = kmer_info['target_var']
    record['target_variant_len'] = len(record['target_variant'])
    record['positive_genotypes'].sort()
    record['positive_genotypes'] = ','.join([str(x) for x in record['positive_genotypes']])
    record['partial_genotypes'].sort()
    record['partial_genotypes'] = ','.join([str(x) for x in record['partial_genotypes']])
    target_variant = record['target_variant']
    is_cds = kmer_info['is_cds']

    if is_cds:
        aa_start = kmer_info['aa_start']+1
        aa_end = kmer_info['aa_end']
        if mutation_type == 'snp':
            record['aa_name'] = "{}{}{}".format(ref_aa, aa_start, target_aa)
        else:
            record['aa_name'] = "{}_{}_{}".format(mutation_type, aa_start, aa_end)

    if mutation_type == 'snp':
        record['mutation_key'] = "{}_{}_{}_{}".format(mutation_type, ref_variant, vStart, alt_variant)
        record['dna_name'] = "{}{}{}".format(ref_variant, vStart, target_variant)
    else:
        record['mutation_key'] = "{}_{}_{}_{}".format(mutation_type, vStart, vEnd, alt_variant)
        record['dna_name'] = "{}_{}_{}_{}".format(mutation_type, vStart, vEnd, target_variant)

    return record

def scheme_format(selected_kmers):
    scheme = {}
    uid = 0
    for mutation_type in selected_kmers:
        for event in selected_kmers[mutation_type]:
            if mutation_type == 'snp':
                vStart = event
                vEnd = event
                for base in selected_kmers[mutation_type][event]['ref']['kmers']:
                    if len(selected_kmers[mutation_type][event]['ref']['kmers'][base]) == 0:
                        continue
                    ref_variant = base
                    break

                for base in selected_kmers[mutation_type][event]['alt']['kmers']:
                    if len(selected_kmers[mutation_type][event]['alt']['kmers'][base]) == 0:
                        continue
                    alt_variant = base
                    for kmer in selected_kmers[mutation_type][event]['ref']['kmers'][ref_variant]:
                        ref_aa = ''
                        if 'ref_var_aa' in selected_kmers[mutation_type][event]['ref']['kmers'][ref_variant][kmer]:
                            ref_aa = selected_kmers[mutation_type][event]['ref']['kmers'][ref_variant][kmer]['ref_var_aa']

                        kmer_info = selected_kmers[mutation_type][event]['ref']['kmers'][ref_variant][kmer]
                        if not 'target_var' in kmer_info:
                            continue
                        scheme[uid] = populate_fields(uid, mutation_type, vStart, vEnd, 'ref', ref_variant, alt_variant,ref_aa,ref_aa,
                                                      kmer_info)
                        uid += 1

                    for kmer in selected_kmers[mutation_type][event]['alt']['kmers'][base]:
                        kmer_info = selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]
                        if not 'target_var' in kmer_info:
                            continue
                        target_aa = ''
                        if 'alt_var_aa' in selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]:
                            target_aa = selected_kmers[mutation_type][event]['alt']['kmers'][base][kmer]['alt_var_aa']
                        scheme[uid] = populate_fields(uid, mutation_type, vStart, vEnd, 'alt', ref_variant, alt_variant,ref_aa,target_aa,
                                                      kmer_info)
                        uid += 1

            else:
                (vStart, vEnd) = event.split(':')
                vStart = int(vStart)
                vEnd = int(vEnd)
                vLen = vEnd - vStart + 1
                if mutation_type == 'del':
                    ref_variant = ''.join(['N'] * vLen)
                    alt_variant = ''.join(['-'] * vLen)
                else:
                    ref_variant = ''.join(['-'] * vLen)
                    alt_variant = ''.join(['N'] * vLen)

                for state in selected_kmers[mutation_type][event]:
                    for kmer in selected_kmers[mutation_type][event][state]['kmers']:
                        kmer_info = selected_kmers[mutation_type][event][state]['kmers'][kmer]
                        if not 'target_var' in kmer_info:
                            continue
                        scheme[uid] = populate_fields(uid, mutation_type, vStart, vEnd, state, ref_variant, alt_variant,'','',
                                                      kmer_info)
                        uid += 1

    return scheme

def get_genotype_kmer_frac(genotype_kCounts, genotype_counts, min_frac):
    rules = {}
    for kmer in genotype_kCounts:
        for genotype in genotype_kCounts[kmer]['counts']:
            if not genotype in genotype_counts:
                continue
            total = genotype_counts[genotype]
            kCount = genotype_kCounts[kmer]['counts'][genotype]
            if total == 0:
                continue
            frac = kCount / total
            if frac > min_frac:
                if not kmer in rules:
                    rules[kmer] = {}
                rules[kmer][genotype] = frac
    return rules

def associate_genotype_rules(selected_kmers, genotype_fracs, min_ref_frac, min_alt_frac):
    min_par_frac = max([1 - min_alt_frac, 1 - min_ref_frac])
    for mutation_type in selected_kmers:
        for event in selected_kmers[mutation_type]:
            if mutation_type == 'snp':
                for state in selected_kmers[mutation_type][event]:
                    thresh = min_ref_frac
                    if state == 'alt':
                        thresh = min_alt_frac
                    for base in selected_kmers[mutation_type][event][state]['kmers']:
                        if len(selected_kmers[mutation_type][event][state]['kmers'][base]) == 0:
                            continue
                        for kmer in selected_kmers[mutation_type][event][state]['kmers'][base]:
                            if not kmer in genotype_fracs:
                                continue
                            gfracs = genotype_fracs[kmer]
                            for genotype in gfracs:
                                frac = gfracs[genotype]
                                if frac >= thresh:
                                    selected_kmers[mutation_type][event][state]['kmers'][base][kmer][
                                        'positive_genotypes'].append(genotype)
                                elif frac >= min_par_frac:
                                    selected_kmers[mutation_type][event][state]['kmers'][base][kmer][
                                        'partial_genotypes'].append(genotype)
            else:
                for state in selected_kmers[mutation_type][event]:
                    thresh = min_ref_frac
                    if state == 'alt':
                        thresh = min_alt_frac

                    for kmer in selected_kmers[mutation_type][event][state]['kmers']:
                        if not kmer in genotype_fracs:
                            continue
                        gfracs = genotype_fracs[kmer]
                        for genotype in gfracs:
                            frac = gfracs[genotype]
                            if frac >= thresh:
                                selected_kmers[mutation_type][event][state]['kmers'][kmer]['positive_genotypes'].append(
                                    genotype)
                            elif frac >= min_par_frac:
                                selected_kmers[mutation_type][event][state]['kmers'][kmer]['partial_genotypes'].append(
                                    genotype)

def write_genotype_reports(out_dir, genotype_counts, genotype_kCounts, scheme):
    mutations_report_file = os.path.join(out_dir, "genotype.mutations.txt")
    kmer_report_file = os.path.join(out_dir, "genotype.kmers.txt")
    inf_kmer_report_file = os.path.join(out_dir, "informative.kmers.txt")

    assoc = {}
    kmer_report = {}
    mutations_report = {}
    for genotype in genotype_counts:
        assoc[genotype] = {
            'total': genotype_counts[genotype],
            'mutations': {'ref': set(), 'alt': set()},
            'kmers': {'ref': {}, 'alt': {}}
        }

        kmer_report[genotype] = {'genotype': genotype,
                                 'num_members': assoc[genotype]['total'],
                                 'num_pos_ref': 0,
                                 'num_pos_alt': 0,
                                 'alt_kmers': [],
                                 'num_uniq': 0,
                                 'uniq_kmers': {}
                                 }

        mutations_report[genotype] = {'genotype': genotype,
                                      'num_members': assoc[genotype]['total'],
                                      'num_pos_ref': 0,
                                      'num_pos_alt': 0,
                                      'alt_mutations': [],
                                      'num_uniq': 0,
                                      'uniq_mutations': []
                                      }

    for uid in scheme:
        state = scheme[uid]['state']
        pos = scheme[uid]['positive_genotypes'].split(',')
        dna_name = scheme[uid]['dna_name']
        kseq = scheme[uid]['kseq']
        is_diag = False
        if len(pos) == 1:
            is_diag = True

        for genotype in pos:
            if not genotype in assoc:
                continue
            assoc[genotype]['mutations'][state].add(dna_name)
            assoc[genotype]['kmers'][state][uid] = kseq
            if state == 'alt':
                mutations_report[genotype]['num_pos_alt'] += 1
                mutations_report[genotype]['alt_mutations'].append(dna_name)
            else:
                mutations_report[genotype]['num_pos_ref'] += 1

            if is_diag:
                mutations_report[genotype]['uniq_mutations'].append(dna_name)
                mutations_report[genotype]['num_uniq'] += 1
                kmer_report[genotype]['uniq_kmers'][uid] = kseq
                kmer_report[genotype]['num_uniq'] += 1

  #  info_kmers = {}
  #  for genotype in assoc:
  #      info_kmers[genotype] = {}
  #      for state in assoc[genotype]['kmers']:
  #          for kmer in assoc[genotype]['kmers'][state]:
  #              info_kmers[genotype][kmer] = calc_kmer_associations(genotype, kmer, genotype_kCounts, genotype_counts, 100)

    #pd.DataFrame.from_dict(info_kmers,orient='index').to_csv(inf_kmer_report_file,sep="\t",header=True, index=False)


    pd.DataFrame.from_dict(kmer_report, orient='index').to_csv(kmer_report_file, sep="\t", header=True, index=False)
    pd.DataFrame.from_dict(mutations_report, orient='index').to_csv(mutations_report_file, sep="\t", header=True, index=False)


    return

def group_kmers_by_start_pos(kmer_results):
    grouped_kmers = {}
    for kIndex in kmer_results:
        aln_start = kmer_results[kIndex]['aStart']
        if not aln_start in grouped_kmers:
            grouped_kmers[aln_start] = {}
        grouped_kmers[aln_start][kIndex] = kmer_results[kIndex]
    return grouped_kmers

def run():
    cmd_args = parse_args()
    logger = init_console_logger(2)
    # input parameters
    input_alignment = cmd_args.input_msa
    input_meta = cmd_args.input_meta
    prefix = cmd_args.prefix
    ref_id = cmd_args.ref_id
    ref_gbk = cmd_args.ref_gbk
    kLen = cmd_args.kmer_len
    out_dir = cmd_args.outdir
    max_ambig = cmd_args.max_ambig
    min_ref_frac = cmd_args.min_ref_frac
    min_alt_frac = cmd_args.min_alt_frac
    max_homo = cmd_args.max_homo
    n_threads = cmd_args.n_threads
    no_plots = cmd_args.no_plots
    min_kmer_count = cmd_args.min_kmer_freq
    min_var_freq = cmd_args.min_var_freq
    batch_size = cmd_args.seq_batch_size

    num_stages = 5
    stage = 0

    logging.info("Parsityper creator v {}".format(__version__))
    if not os.path.isdir(out_dir):
        logging.info("Creating analysis directory {}".format(out_dir))
        os.mkdir(out_dir, 0o755)
    if max_homo is None:
        max_homo = kLen


    logging.info("Performing parsityper creator analysis")


    if min_alt_frac < min_ref_frac:
        min_frac = min_alt_frac
    else:
        min_frac = min_ref_frac

    total_sys_memory = psutil.virtual_memory().total
    total_sys_threads = psutil.cpu_count()
    if n_threads > total_sys_threads:
        n_threads = total_sys_threads

    # read metadata
    logging.info("Reading genotype assignments from {}".format(input_meta))
    metadata_df = read_tsv(input_meta)
    logging.info("Found {} lines in {}".format(len(metadata_df), input_meta))
    metadata_df['sample_id'] = metadata_df['sample_id'].astype(str)
    metadata_df['genotype'] = metadata_df['genotype'].astype(str)
    genotype_mapping = get_genotype_mapping(metadata_df)

    if batch_size is None or batch_size == 0:
        batch_size = int(len(genotype_mapping) / n_threads)

    # Get the Gene features from the reference sequenc if specified
    perform_annotation = False
    ref_features = {}

    if ref_gbk is not None:
        perform_annotation = True
        if os.path.isfile(ref_gbk):
            ref_features = parse_reference_sequence(ref_gbk)
        else:
            logging.error("Specified ref_gbk does not exist: {}".format(ref_gbk))

    # process msa
    stime = time.time()
    fasta_dir = os.path.join(out_dir, "_stage-{}".format(stage))
    if not os.path.isdir(fasta_dir):
        logging.info("Creating analysis directory for interim fastas {}".format(fasta_dir))
        os.mkdir(fasta_dir, 0o755)

    logging.info("stage-{}: Processing sequences from MSA {}".format(stage,input_alignment))
    msa_info = preprocess_seqs(input_alignment, ref_id, list(genotype_mapping.keys()), batch_size, fasta_dir, 'stage-0',
                               iupac_replacement, min_var_freq)
    if not msa_info['is_align_ok']:
        logging.error("Input alignment has issues, please correct and try again")
        sys.exit()
    align_len = msa_info['align_len']
    num_sequences = len(msa_info['seq_ids'])
    max_kmer_count = num_sequences


    stage+=1

    logging.info("stage-{}: Filtering samples from metadata which do not have a sequence in MSA".format(stage))
    genotype_counts = {}
    filter_samples = {}
    for sample_id in genotype_mapping:
        genotype = genotype_mapping[sample_id]
        if sample_id in msa_info['seq_ids']:
            filter_samples[sample_id] = genotype
            if not genotype in genotype_counts:
                genotype_counts[genotype] = 0
            genotype_counts[genotype] += 1
        else:
            logging.warning("stage-{}: sample {} does not have a sequence in MSA".format(stage,sample_id))
    genotype_mapping = filter_samples
    num_samples = len(genotype_mapping)
    del (filter_samples)

    # kmer counting
    init_jellyfish_mem = int(num_samples * batch_size / 1000000)
    if init_jellyfish_mem == 0:
        init_jellyfish_mem = 1
    logging.info("stage-{}: Initial jellyfish cache size set to {}M".format(stage,init_jellyfish_mem))
    logging.info("stage-{}: Perfoming k-mer counting using {} threads".format(stage, n_threads))
    perform_kmer_counting(msa_info['unalign_files'], kLen, "{}M".format(init_jellyfish_mem), n_threads)
    logging.info("stage-{}: Combining {} jellyfish files".format(stage, len(msa_info['unalign_files'])))
    seqKmers = filter_kmers(combine_jellyfish_results(msa_info['unalign_files']), min_kmer_count, max_kmer_count,
                            max_ambig, max_homo)

    logging.info("stage-{}: Perfoming k-mer searching using {} threads".format(stage, n_threads))
    aho_results = map_kmers(seqKmers, msa_info['subset_files'], n_threads)
    add_ref_kmer_info(aho_results, seqKmers, msa_info['ref_aln'], kLen)

    logging.info("stage-{}: Grouping k-mers by start position".format(stage))
    kmer_alignment_positions = group_kmers_by_start_pos(aho_results['kmer'])
    logging.info("stage-{}: Filtering singleton k-mers".format(stage))
    filt = {}
    for pos in kmer_alignment_positions:
        if len(kmer_alignment_positions[pos]) == 1:
            continue
        filt[pos] = kmer_alignment_positions[pos]
    kmer_alignment_positions = filt

    logging.info("stage-{}: Selecting kmers".format(stage))
    selected_kmers = {
        'snp': select_snp_kmers(msa_info['variants']['snp'], kLen, kmer_alignment_positions),
        'del': select_indel_kmers(msa_info['variants']['del'], kLen, kmer_alignment_positions),
        'ins': select_indel_kmers(msa_info['variants']['ins'], kLen, kmer_alignment_positions)}

    kmer_indicies = getSelectedKmerIndicies(selected_kmers)
    logging.info("stage-{}: {} kmers selected".format(stage,len(kmer_indicies)))
    logging.info("stage-{}: associating {} kmers with {} genotypes".format(stage,len(kmer_indicies),len(genotype_counts)))
    genotype_kCounts = get_kmer_genotype_counts(kmer_indicies, aho_results['samples'], genotype_mapping)
    del (seqKmers)

    stage+=1

    logging.info("stage-{}: Calculating kmer entropy by genotype".format(stage))
    kmer_entropies = calc_kmer_entropy(genotype_kCounts)


    # Add in entropy and count information
    for mutation_type in selected_kmers:
        for pos in selected_kmers[mutation_type]:
            for state in selected_kmers[mutation_type][pos]:
                if mutation_type == 'snp':
                    for base in selected_kmers[mutation_type][pos][state]['kmers']:
                        for kmer_id in selected_kmers[mutation_type][pos][state]['kmers'][base]:
                            selected_kmers[mutation_type][pos][state]['kmers'][base][kmer_id]['entropy'] = \
                            kmer_entropies[kmer_id]
                            selected_kmers[mutation_type][pos][state]['kmers'][base][kmer_id]['total_count'] = \
                            genotype_kCounts[kmer_id]['total']
                else:
                    for kmer_id in selected_kmers[mutation_type][pos][state]['kmers']:
                        selected_kmers[mutation_type][pos][state]['kmers'][kmer_id]['entropy'] = kmer_entropies[kmer_id]
                        selected_kmers[mutation_type][pos][state]['kmers'][kmer_id]['total_count'] = \
                        genotype_kCounts[kmer_id]['total']

    stage+=1

    logging.info("stage-{}: Calculating kmer positivity rate".format(stage))
    genotype_kmer_fracs = get_genotype_kmer_frac(genotype_kCounts, genotype_counts, min_frac)

    logging.info("stage-{}: Adding genotyping rules".format(stage))
    associate_genotype_rules(selected_kmers, genotype_kmer_fracs, min_ref_frac, min_alt_frac)
    add_gene_inference(selected_kmers, msa_info['ref_aln'], ref_id, ref_features, trans_table=1)


    stage+=1

    logging.info("stage-{}: Formatting scheme".format(stage))
    scheme = scheme_format(selected_kmers)

    logging.info("Writting scheme file")
    select_kmer_file = os.path.join(out_dir, "{}-scheme.txt".format(prefix))
    pd.DataFrame.from_dict(scheme, orient='index').to_csv(select_kmer_file, sep="\t", header=True, index=False)

    logging.info("stage-{}: Determining information content of each scheme kmer".format(stage))
    write_genotype_reports(out_dir, genotype_counts, genotype_kCounts, scheme)

    logging.info("Cleaning up interim files")
    for i in range(0,stage):
        tmp =os.path.join(out_dir, "stage-{}.pickle".format(stage))
        if os.path.isfile(tmp):
            os.remove(tmp)
        tmp = os.path.join(out_dir, "_stage-{}".format(stage))
        if os.path.isdir(tmp):
            os.rmdir(tmp)
    logging.info("Run complete")
