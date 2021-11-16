#!/usr/bin/python
import time
from argparse import (ArgumentParser)
import logging, os, sys, operator
from collections import Counter
import pandas as pd
from parsityper.helpers import init_console_logger, read_tsv, parse_reference_sequence, calc_consensus, generate_consensus_seq,\
find_snp_positions,  find_internal_gaps, count_kmers
from parsityper.helpers import  read_fasta, get_aa_delta, generate_non_gap_position_lookup
from parsityper.visualizations import dendrogram_visualization
from parsityper.scheme import SCHEME_HEADER, parseScheme, constructSchemeLookups, detectAmbigGenotypes
from parsityper.kmerSearch.kmerSearch import init_automaton_dict, perform_kmerSearch_fasta,process_kmer_results
from parsityper.ext_tools.mash import mash_sample_comparisons
from multiprocessing import Pool
from parsityper.ext_tools.jellyfish import run_jellyfish_count,parse_jellyfish_counts
from parsityper.kmerSearch.kmerSearch import  revcomp


def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Kmer scheme generator')
    parser.add_argument('--input_msa', type=str, required=True,
                        help='MSA in fasta format')
    parser.add_argument('--input_meta', type=str, required=True,
                        help='TSV file of sample_id,genotype')
    parser.add_argument('--ref_id', type=str, required=True,
                        help='sample_id for reference sequence to use in MSA')
    parser.add_argument('--ref_gbk', type=str, required=True,
                        help='GenBank file for reference sequences')
    parser.add_argument('--outdir', type=str, required=True,
                        help='output directory')
    parser.add_argument('--prefix', type=str, required=False,
                        help='output file prefix',default='parsityper')
    parser.add_argument('--min_len', type=int, required=False,
                        help='Absolute minimum length of acceptable k-mer',default=18)
    parser.add_argument('--max_len', type=int, required=False,
                        help='Absolute minimum length of acceptable k-mer',default=18)
    parser.add_argument('--max_ambig', type=int, required=False,
                        help='Absolute maximum of degenerate bases allowed in a k-mer',default=10)
    parser.add_argument('--max_states', type=int, required=False,
                        help='Absolute maximum of states allowed per kmer',default=256)
    parser.add_argument('--min_complexity', type=int, required=False,
                        help='Absolute maximum of dimer composition',default=0.6)
    parser.add_argument('--max_missing', type=float, required=False,
                        help='Absolute maximum percentage of sequences allowed to be missing kmer',default=0.25)
    parser.add_argument('--n_threads', type=int, required=False,
                        help='Number of threads to use',default=1)
    parser.add_argument('-q', '--impute', required=False, help='impute ambiguous bases',
                        action='store_true')
    parser.add_argument('--iFrac', type=float, required=False,
                        help='fraction of bases needed for imputing ambiguous bases',default=0.9)
    return parser.parse_args()


def get_non_gap_position(ref_non_gap_lookup,pos):

    non_gap_position = ref_non_gap_lookup[pos]
    while non_gap_position == -1:
        pos -= 1
        non_gap_position = ref_non_gap_lookup[pos]
    return non_gap_position

def runKmerCounting(input_seqs,out_dir,kLen,n_threads=1):
    kCount_files = []
    seqFiles = []
    for seq_id in input_seqs:
        seq_file = os.path.join(out_dir, "{}.fasta".format(seq_id))
        out_file = os.path.join(out_dir, "{}.mers".format(seq_id))
        kCount_files.append(out_file)
        seqFiles.append(seq_file)
        fh = open(seq_file, 'w')
        fh.write(">{}\n{}\n".format(seq_id, input_seqs[seq_id]))
        fh.close()
    if n_threads > 1:
        pool = Pool(processes=n_threads)
        res = []
        for i in range(0,len(seqFiles)):
           res.append(pool.apply_async(run_jellyfish_count,(seqFiles[i],kCount_files[i],kLen,n_threads)))
        pool.close()
        pool.join()
        for i in range(0, len(res)):
            res[i].get()
    else:
       for i in range(0,len(seqFiles)):
           run_jellyfish_count(seqFiles[i],kCount_files[i],kLen,n_threads)

    return kCount_files

def getKmerCounts(input_seqs,out_dir,kLen,max_count=1,max_ambig=0,n_threads=1):
    stime = time.time()
    kCount_files = runKmerCounting(input_seqs,out_dir,kLen,n_threads)
    stime = time.time()
    kMers = {}
    for file in kCount_files:
        sample_id = os.path.basename(file).replace('.mers','')
        kmer_df = parse_jellyfish_counts(file)
        kmer_df = kmer_df[kmer_df['count'] <= max_count]
        for row in kmer_df.itertuples():
            kmer = row.kmer
            bCounts = Counter(kmer)
            if 'N' in bCounts:
                if bCounts['N'] > max_ambig:
                    continue
            if not kmer in kMers:
                kMers[kmer] = []
            kMers[kmer].append(sample_id)
    return kMers

def generate_gap_lookup(aln,seq):
    aLen = len(aln)
    sLen = len(seq)
    seq_lookup = []
    p = 0
    for i in range(0,sLen):
        padding = 0
        for k in range(p,aLen):
            aBase = aln[k]
            sBase = seq[i]
            if aBase == '-':
                padding+=1
                continue
            if aBase == sBase:
                break
        p = k
        seq_lookup.append(p)
    return seq_lookup

def getCandidateKmers(variant_positions,input_alignment,unaligned,unal_to_aln_coords,ref_id,path,kLen,max_count=1,max_ambig=0,n_threads=1):
    raw_kmers = getKmerCounts(unaligned, path, kLen, max_count, max_ambig,
                              n_threads)
    aln_len = len(input_alignment[ref_id])
    num_samples = len(input_alignment)
    filtered_kmers = {}

    for kmer in raw_kmers:
        kLen = len(kmer)
        seq_ids = raw_kmers[kmer]
        if len(seq_ids) == num_samples:
            continue
        seq_id = seq_ids[0]
        kStart_unaln = unaligned[seq_id].find(kmer)

        kEnd_unaln = kStart_unaln + kLen
        kStart_aln = unal_to_aln_coords[seq_id][kStart_unaln]
        # revcomp kmer sequences to correct strand
        if kStart_unaln == -1:
            kmer = revcomp(kmer)
            kStart_unaln = unaligned[seq_id].find(kmer)
            kEnd_unaln = kStart_unaln + kLen
            kStart_aln = unal_to_aln_coords[seq_id][kStart_unaln]

        count_bases = 1
        for i in range(kStart_aln, aln_len):
            if input_alignment[seq_id][i] == '-':
                continue
            count_bases += 1
            if count_bases == kLen:
                break
        kEnd_aln = i
        overlapping_variants = []
        for (vStart, vEnd) in variant_positions:
            if vStart >= kStart_aln and vStart <= kEnd_aln:
                overlapping_variants.append([vStart, vEnd])
            elif vEnd >= kStart_aln and vEnd <= kEnd_aln:
                overlapping_variants.append([vStart, vEnd])

        # disregard kmers which are not involved in a variant
        if len(overlapping_variants) == 0:
            continue

        filtered_kmers[kmer] = {'kStart_unalign': kStart_unaln,
                                'kEnd_unalign': kEnd_unaln,
                                'kStart_align': kStart_aln,
                                'kEnd_align': kEnd_aln,
                                'affect_variants': overlapping_variants,
                                'num_seqs': len(seq_ids),
                                'seq_ids': seq_ids}
    return filtered_kmers

def build_mutation_lookup(variant_positions,ref_id,input_alignment,kmer_len):
    # init ref kmer data structure
    ref_seq = input_alignment[ref_id]
    ref_len = len(ref_seq)
    ref_seq_nmasked = ref_seq.replace('A', 'N').replace('T', 'N').replace('C','N').replace('G', 'N')
    num_variants = len(variant_positions)
    mutations = {}
    for i in range(0, num_variants):
        vStart = variant_positions[i][0]
        vEnd = variant_positions[i][1]
        kStart = vStart - kmer_len
        kEnd = vEnd + 1
        if kStart < 0:
            kStart = 0
        if kEnd >= ref_len:
            kEnd = ref_len - 1

        variant = ref_seq_nmasked[vStart:vEnd +1]
        if vStart - vEnd == 0 and variant == 'N':
            mutation_type = 'snp'
            mutation_key = "snp{}".format(vStart)
            variant = ref_seq[vStart]
        else:
            if '-' not in variant:
                mutation_type = 'del'
                mutation_key = "del{}_{}".format(vStart,vEnd)
            else:
                mutation_type = 'ins'
                mutation_key = "ins{}_{}".format(vStart, vEnd)

        if not mutation_key in mutations:
            mutations[mutation_key] = {'mutation_type':mutation_type,'vStart':vStart,
                                       'vEnd':vEnd,'vLen':(vEnd-vStart)+1,'kStart':kStart,'kEnd':kEnd,'ref_variant':variant}

    return mutations

def get_alt_variants(input_alignment,mutations,align_len):
    variants = {}
    #process each sequence to create a list of the variants
    for seq_id in input_alignment:
        for mutation_key in mutations:
            vStart = mutations[mutation_key]['vStart']
            vEnd = mutations[mutation_key]['vEnd']
            mutation_type = mutations[mutation_key]['mutation_type']
            ref_variant = mutations[mutation_key]['ref_variant']
            variant = input_alignment[seq_id][vStart:vEnd + 1]

            if mutation_type == 'snp':
                if variant == 'N':
                    continue
                if '-' == variant:
                    if vStart > 0:
                        s = vStart - 1
                    if vEnd < len(input_alignment[seq_id]) -2:
                        e = vEnd +1
                    if input_alignment[seq_id][s] == '-' and input_alignment[seq_id][e] == '-':
                        continue
                if variant == ref_variant:
                    state = 'ref'
                else:
                    state = 'alt'
            else:
                # Anchor variant to an N base around gaps
                # this handels multiple indel substrings
                iStart = vStart
                while input_alignment[seq_id][iStart] == '-' and iStart > 0:
                    iStart -= 1
                iEnd = vEnd
                while input_alignment[seq_id][iEnd] == '-' and iEnd < align_len -1:
                    iEnd += 1
                if mutation_type == 'ins':
                    if abs(iStart - vStart) == 0 and (iEnd - vEnd) == 0:
                        state = 'alt'
                    else:
                        state = 'ref'
                else:
                    if abs(iStart - vStart) == 1 and (iEnd - vEnd) == 1:
                        state = 'alt'
                    else:
                        state = 'ref'
            #if state == 'alt':
            #    print("{}\t{}\t{}\t{}\t{}\t{}".format(seq_id,state,vStart,vEnd,variant,ref_variant))
            if state == 'ref':
                continue
            if not mutation_key in variants:
                variants[mutation_key] = {}
            if not variant in variants[mutation_key]:
                variants[mutation_key][variant] = {
                    'mutation_type': mutation_type,
                    'vStart': vStart,
                    'vEnd': vEnd,
                    'vLen': (vEnd - vStart) + 1,
                    'state':state,
                    'target_variant': variant,
                    'seq_ids':[]
                }
            variants[mutation_key][variant]['seq_ids'].append(seq_id)

    return variants

def select_overlapping_kmers(raw_kmers,vStart,vEnd):
    filt_kmers = {}
    for kmer in raw_kmers:
        kStart_aln = raw_kmers[kmer]['kStart_align']
        kEnd_aln = raw_kmers[kmer]['kEnd_align']
        is_found = False
        if vStart >= kStart_aln and vStart <= kEnd_aln:
            is_found=True
        elif vEnd >= kStart_aln and vEnd <= kEnd_aln:
            is_found=True
        if is_found:
            filt_kmers[kmer] = raw_kmers[kmer]
    return filt_kmers

def create_entry(variant_key,dna_name,mutation_type,groups,vStart,vEnd,unalign_vstart,unalign_vend,variant_seq,ref_seq,ref_non_gap_lookup,state):
    kmer_entries = []
    for kmer in groups:
        kStart = groups[kmer]['kStart_align']
        kEnd = groups[kmer]['kEnd_align']
        kMembers = groups[kmer]['seq_ids']
        unalign_kstart = get_non_gap_position(ref_non_gap_lookup, kStart)
        unalign_kend = get_non_gap_position(ref_non_gap_lookup, kEnd)
        kmer_entries.append(
            {
                'mutation_key': variant_key,
                'dna_name': dna_name,
                'align_variant_start': vStart,
                'align_variant_end': vEnd,
                'unalign_variant_start': unalign_vstart,
                'unalign_variant_end': unalign_vend,
                'align_kmer_start': kStart,
                'align_kmer_end': kEnd,
                'unalign_kmer_start': unalign_kstart,
                'unalign_kmer_end': unalign_kend,
                'target_variant': variant_seq,
                'target_variant_len': len(variant_seq),
                'ref_variant': ref_seq,
                'ref_variant_len': len(ref_seq),
                'mutation_type': mutation_type,
                'state': state,
                'align_kseq': kmer,
                'unalign_kseq': kmer,
                'unalign_klen': len(kmer),
                'seq_ids': kMembers,
                'positive_genotypes': [],
                'partial_genotypes': [],
                'negative_genotypes': [],
                'is_kmer_found': True,
                'is_kmer_length_ok': True,
                'is_kmer_unique': True,
                'is_kmer_complexity_ok': True,
                'is_valid': True
            })
    return kmer_entries

def process_variant(kmers,mutation_type,vStart,vEnd,unalign_vstart,unalign_vend,ref_non_gap_lookup,ref_variant_seq,variant_events,seq_ids,seq_kmers,min_kmer_len, max_kmer_len):
    kmer_scheme = {}
    #get kmers exclusive to each variant
    for variant in variant_events:
        in_seq_ids = variant_events[variant]['seq_ids']
        out_seq_ids = list(set(seq_ids) - set(in_seq_ids))
        valid_alt_kmers = {}
        valid_ref_kmers = {}
        for kmer in kmers:
            kmer_seq_ids = list(set(kmers[kmer]['seq_ids']))
            if len(kmer_seq_ids) == 0:
                continue
            in_overlap = len(set(in_seq_ids) & set(kmer_seq_ids))
            if in_overlap == len(kmer_seq_ids):
                valid_alt_kmers[kmer] = len(kmer_seq_ids)
            elif in_overlap == 0:
                valid_ref_kmers[kmer] = len(kmer_seq_ids)

        if len(valid_alt_kmers) == 0:
            continue

        valid_alt_kmers = sorted(valid_alt_kmers.items(), key=lambda x: x[1], reverse=True)
        valid_ref_kmers = sorted(valid_ref_kmers.items(), key=lambda x: x[1], reverse=True)

        selected_alt_kmers = {}
        assigned_alt_samples = []
        num_in_samples = len(in_seq_ids)
        for (kmer,count) in valid_alt_kmers:
            kmer_seq_ids = list(set(kmers[kmer]['seq_ids']))
            unassigned_samples = list(set(kmer_seq_ids) - set(assigned_alt_samples))
            if len(unassigned_samples) == 0:
                continue
            assigned_alt_samples.extend(unassigned_samples)
            selected_alt_kmers[kmer] = kmers[kmer]
            if len(assigned_alt_samples) == num_in_samples:
                break

        selected_ref_kmers = {}
        assigned_ref_samples = []
        num_out_samples = len(out_seq_ids)
        for (kmer,count) in valid_ref_kmers:
            kmer_seq_ids = list(set(kmers[kmer]['seq_ids']))
            unassigned_samples = list(set(kmer_seq_ids) - set(assigned_ref_samples))
            if len(unassigned_samples) == 0:
                continue
            assigned_ref_samples.extend(unassigned_samples)
            selected_ref_kmers[kmer] = kmers[kmer]
            if len(assigned_ref_samples) == num_out_samples:
                break

        if mutation_type == 'snp':
            variant_key = "{}_{}_{}_{}".format(mutation_type, ref_variant_seq, vStart, variant)
            dna_name = "{}{}{}".format(ref_variant_seq, vStart, variant)
        else:
            variant_key = "{}_{}_{}".format(mutation_type, vStart, vEnd)
            dna_name = variant_key
        kmer_scheme[variant_key] = {'ref': [], 'alt': []}
        kmer_scheme[variant_key]['alt'] = create_entry(variant_key, dna_name, mutation_type, selected_alt_kmers , vStart, vEnd,
                                                       unalign_vstart, unalign_vend,
                                                       variant, ref_variant_seq, ref_non_gap_lookup, 'alt')
        kmer_scheme[variant_key]['ref'] = create_entry(variant_key, dna_name, mutation_type, selected_ref_kmers , vStart,
                                                       vEnd, unalign_vstart, unalign_vend,
                                                       variant, ref_variant_seq, ref_non_gap_lookup, 'ref')
    return kmer_scheme

def calc_kmer_complexity(kSeq):
    mers = count_kmers(kSeq, K=2)
    num_kmers = sum(mers.values())
    mer_perc = []
    for m in mers:
        mer_perc.append(mers[m] / num_kmers)
    if len(mer_perc) > 0:
        m = max(mer_perc)
    else:
        m = 1
    return m

def calc_complexity_scheme(scheme):
    for mutation in scheme:
        states = ['ref','alt']
        for state in states:
            for entry in scheme[mutation][state]:
                kSeq = entry['unalign_kseq']
                entry['complexity'] = calc_kmer_complexity(kSeq)
    return scheme

def construct_scheme(variant_positions, ref_id, input_alignment, jellyfish_path,min_kmer_len,max_kmer_len,n_threads=1):
    # initialize analysis directory
    if not os.path.isdir(jellyfish_path):
        logging.info("Creating analysis results directory {}".format(jellyfish_path))
        os.mkdir(jellyfish_path, 0o755)
    else:
        logging.info("Results directory {} already exits, will overwrite any results files here".format(jellyfish_path))

    seq_ids = list(input_alignment.keys())
    ref_non_gap_lookup = generate_non_gap_position_lookup(input_alignment[ref_id])
    unaligned_seqs = {}
    unal_to_aln_coords = {}
    logging.info("Building unaligned lookup")
    for seq_id in input_alignment:
        unaligned_seqs[seq_id] = input_alignment[seq_id].replace('-','')
        unal_to_aln_coords[seq_id] = generate_gap_lookup(input_alignment[seq_id], unaligned_seqs[seq_id])



    align_len = len(input_alignment[ref_id])
    logging.info("Building mutation event lookup")
    mutations = build_mutation_lookup(variant_positions, ref_id, input_alignment, min_kmer_len)
    logging.info("Cataloguing variants")
    variant_events = get_alt_variants(input_alignment,mutations,align_len)
    positions = []

    for mutation_key in variant_events:
        if not mutation_key in variant_events:
            logging.info("Mutation {} did not have a valid alt state..skipping".format(mutation_key))
            continue
        for variant in variant_events[mutation_key]:
            vStart = variant_events[mutation_key][variant]['vStart']
            vEnd = variant_events[mutation_key][variant]['vEnd']
            positions.append([vStart,vEnd])

    if n_threads > 1:
        pool = Pool(processes=n_threads)

    logging.info("Begining k-mer extraction for {} variant positions".format(len(positions)))
    stime = time.time()

    seq_kmers = getCandidateKmers(variant_positions,input_alignment,unaligned_seqs,unal_to_aln_coords,ref_id,jellyfish_path,max_kmer_len,max_count=1,max_ambig=0,n_threads=n_threads)
    logging.info("time for kmer extraction {}".format(time.time() - stime))

    stime = time.time()
    res = []
    kmer_scheme = {}
    count = 0
    logging.info("Deconvoluting complex mutations")

    #triage single bp deletions
    mutation_keys = list(mutations.keys())
    for mutation_key in mutation_keys:
        mutation_type = mutations[mutation_key]['mutation_type']
        if mutation_type != 'snp':
            continue
        vStart = mutations[mutation_key]['vStart']
        vEnd = mutations[mutation_key]['vEnd']
        kStart = mutations[mutation_key]['kStart']
        kEnd = mutations[mutation_key]['kEnd']
        variant = mutations[mutation_key]['ref_variant']
        vlen = mutations[mutation_key]['vLen']
        if '-' in variant_events[mutation_key]:
            m = {}
            for key in variant_events[mutation_key]['-']:
                m[key] = variant_events[mutation_key]['-'][key]
            m['mutation_type'] = 'del'
            m['target_variant'] = '-'
            key = "del{}_{}".format(vStart,vEnd)
            variant_events[key] = {}
            variant_events[key]['-'] = m
            mutations[key] = {'mutation_type':'del','vStart':vStart,
                                       'vEnd':vEnd,'vLen':vlen,'kStart':kStart,'kEnd':kEnd,'ref_variant':variant}
            del(variant_events[mutation_key]['-'])
            if len(variant_events[mutation_key]) == 0:
                del (mutations[mutation_key])

    logging.info("Optimizing k-mer selection for {} variant positions using ranges {}-{}".format(len(positions),min_kmer_len,max_kmer_len))
    for mutation_key in mutations:
        count+=1
        if not mutation_key in variant_events:
            continue
        stime = time.time()
        mutation_type = mutations[mutation_key]['mutation_type']
        vStart = mutations[mutation_key]['vStart']
        vEnd = mutations[mutation_key]['vEnd']
        ovKmers = select_overlapping_kmers(seq_kmers, vStart, vEnd)
        if len(ovKmers) == 0:
            continue

        unalign_vstart =  get_non_gap_position(ref_non_gap_lookup, vStart)
        unalign_vend =  get_non_gap_position(ref_non_gap_lookup, vEnd)
        ref_variant_seq = input_alignment[ref_id][vStart:vEnd +1]
        ovKmers = select_overlapping_kmers(seq_kmers,vStart,vEnd)

        if n_threads == 1:
            kmer_scheme.update(process_variant(ovKmers,mutation_type,vStart,vEnd,unalign_vstart, unalign_vend,ref_non_gap_lookup,ref_variant_seq,variant_events[mutation_key],seq_ids,ovKmers,min_kmer_len, max_kmer_len))
        else:
            res.append(pool.apply_async(process_variant, (ovKmers,mutation_type,vStart,vEnd,unalign_vstart, unalign_vend,ref_non_gap_lookup,ref_variant_seq,variant_events[mutation_key],seq_ids,ovKmers,min_kmer_len, max_kmer_len)))

    if n_threads > 1:
        pool.close()
        pool.join()
    for i in range(0, len(res)):
        kmer_scheme.update(res[i].get())
    logging.info("time for kmer selection {}".format(time.time() - stime))
    return kmer_scheme

def add_gene_inference(scheme,ref_name,reference_info,input_alignment):
    for mutation_key in scheme:
        for state in scheme[mutation_key]:
            for record in scheme[mutation_key][state]:
                mutation_type = record['mutation_type']
                align_start = record['align_variant_start']
                align_end= record['align_variant_end']
                non_gap_start = record['unalign_variant_start']
                non_gap_end= record['unalign_variant_end']
                seq_id = record['seq_ids'][0]

                #Fix Variant information
                if mutation_type == 'del':
                    variant = ''.join(['-']*(1+align_end-align_start))
                elif mutation_type == 'ins':
                    variant = input_alignment[seq_id][align_start:align_end+1].replace('-','')
                else:
                    variant = record['target_variant']

                record['target_variant'] = variant
                record['target_variant_len'] = len(variant)
                aa_info = get_aa_delta(non_gap_start, non_gap_end, variant, mutation_type,reference_info,ref_name,trans_table=1)
                for key in aa_info:
                    value = aa_info[key]
                    record[key] = value

    return scheme

def qa_scheme(scheme,input_alignment,ref_id,genotype_mapping,min_len, max_len,min_complexity,n_threads=1):
    scheme = calc_complexity_scheme(scheme)
    unaligned = {}
    for seq_id in input_alignment:
        unaligned[seq_id] = input_alignment[seq_id].replace("-",'')
    scheme = multiplicity_check(scheme, unaligned,genotype_mapping)
    problem_mutations = []
    for mutation_key in scheme:
        for state in scheme[mutation_key]:
            for row in scheme[mutation_key][state]:
                klen = len(row['unalign_kseq'])
                if klen > max_len or klen < min_len:
                    row['is_kmer_length_ok'] = False
                    row['is_valid'] = False
                if row['complexity'] >= 1 - min_complexity:
                    row['is_kmer_complexity_ok'] = False
                    row['is_valid'] = False
                if not row['is_valid']:
                    problem_mutations.append(mutation_key)
    return scheme

def format_scheme_human_readable(scheme,ref_id,input_alignment):
    seq = input_alignment[ref_id]
    for key in scheme:
        for state in scheme[key]:
            for row in scheme[key][state]:
                ref_variant = seq[row['align_variant_start']:row['align_variant_end'] + 1].replace('-', '')
                #change from base 0 to base 1 for all coordinates
                row['align_variant_start']+=1
                row['align_variant_end'] += 1
                row['unalign_variant_start'] += 1
                row['unalign_variant_end'] += 1
                row['align_kmer_start']+=1
                row['align_kmer_end'] += 1
                row['unalign_kmer_start']+=1
                row['unalign_kmer_end'] += 1
                row['gene_start'] += 1
                row['gene_end'] += 1
                row['cds_start'] += 1
                row['cds_end'] += 1
                row['aa_start'] += 1
                row['aa_end'] += 1

                target_variant = row['target_variant']
                mutation_type = row['mutation_type']
                is_cds = row['is_cds']
                aa_name = 'N/A'
                ref_state = row['ref_state']
                alt_state = row['alt_state']
                aa_start = row['aa_start']
                aa_end = row['aa_end']
                if mutation_type == 'snp':
                    mutation_key = "{}_{}".format(mutation_type,row['unalign_variant_start'])
                    if state == 'alt':
                        dna_name = "{}{}{}".format(ref_variant,row['unalign_variant_start'],target_variant)
                    else:
                        dna_name = "{}{}{}".format(ref_variant, row['unalign_variant_start'], ref_variant)
                    if is_cds:
                        aa_name = "{}{}{}".format(ref_state,aa_start,ref_state)
                        if state == 'alt':
                            aa_name = "{}{}{}".format(ref_state,aa_start,alt_state)
                elif mutation_type == 'del':
                    mutation_key = "{}_{}_{}".format(mutation_type,row['unalign_variant_start'],row['unalign_variant_end'])
                    dna_name = "{}_{}_{}_{}".format(mutation_type,row['unalign_variant_start'],row['unalign_variant_end'],ref_variant )
                    if state == 'alt':
                        dna_name = "{}_{}_{}".format(mutation_type, row['unalign_variant_start'],
                                                        row['unalign_variant_end'])
                    if is_cds:
                        aa_name = "{}_{}_{}_{}".format(mutation_type, aa_start, aa_end, ref_state)
                        if state == 'alt':
                            aa_name = "{}_{}_{}".format(mutation_type, aa_start, aa_end)
                else:
                    mutation_key = "{}_{}_{}".format(mutation_type, row['unalign_variant_start'],row['unalign_variant_end'])
                    if state == 'alt':
                        dna_name = "{}_{}_{}_{}".format(mutation_type,row['unalign_variant_start'],row['unalign_variant_end'],target_variant )
                    else:
                        dna_name = "{}_{}_{}".format(mutation_type, row['unalign_variant_start'],
                                                        row['unalign_variant_end'])
                    if is_cds:
                        if state == 'alt':
                            aa_name = "{}_{}_{}_{}".format(mutation_type, aa_start, aa_end, alt_state)
                        else:
                            aa_name = "{}_{}_{}".format(mutation_type, aa_start, aa_end)
                row['mutation_key'] = mutation_key
                row['dna_name'] = dna_name
                row['aa_name'] = aa_name

    return  scheme

def identify_shared_kmers(genotype_mapping,scheme,min_thresh=0.5,max_thresh=1):
    '''
    :param genotype_mapping:
    :type genotype_mapping:
    :param kmer_profile:
    :type kmer_profile:
    :param min_thresh:
    :type min_thresh:
    :param max_thresh:
    :type max_thresh:
    :return:
    :rtype:
    '''

    all_seq_ids = list(genotype_mapping.keys())

    #Get counts of each genotype
    genotype_sample_counts = {}
    for sample_id in genotype_mapping:
        genotype = genotype_mapping[sample_id]
        if not genotype in genotype_sample_counts:
            genotype_sample_counts[genotype] = 0
        genotype_sample_counts[genotype]+=1
    genotype_shared_kmers = {}
    for genotype in genotype_sample_counts:
        genotype_shared_kmers[genotype] = []


    #Process each kmer and determine if it meets the thresholds

    for mutation_key in scheme:
        for state in scheme[mutation_key]:
            for kmer_entry in scheme[mutation_key][state]:
                uid =  kmer_entry['key']
                seq_ids = kmer_entry['seq_ids']
                genotypes = []
                for sample_id in seq_ids:
                    genotype = genotype_mapping[sample_id]
                    genotypes.append(genotype)
                gCounts = Counter(genotypes)
                for genotype in gCounts:
                    count = gCounts[genotype]
                    perc = count / genotype_sample_counts[genotype]
                    if perc >= min_thresh and perc <= max_thresh:
                        genotype_shared_kmers[genotype].append(uid)
    return genotype_shared_kmers

def identify_shared_mutations(genotype_mapping,scheme,min_thresh=0.5,max_thresh=1):
    '''
    :param genotype_mapping:
    :type genotype_mapping:
    :param kmer_profile:
    :type kmer_profile:
    :param min_thresh:
    :type min_thresh:
    :param max_thresh:
    :type max_thresh:
    :return:
    :rtype:
    '''
    #Get counts of each genotype
    genotype_sample_counts = {}
    shared_mutations = {}
    for sample_id in genotype_mapping:
        genotype = genotype_mapping[sample_id]
        if not genotype in genotype_sample_counts:
            genotype_sample_counts[genotype] = 0
            shared_mutations[genotype] = []
        genotype_sample_counts[genotype]+=1

    #get counts of genotypes by mutation
    mutation_geno_counts = {}
    for mutation_key in scheme:
        for state in scheme[mutation_key]:
            for row in scheme[mutation_key][state]:
                dna_name = row['dna_name']
                if not dna_name in mutation_geno_counts:
                    mutation_geno_counts[dna_name] = {}
                seq_ids = row['seq_ids']
                for sid in seq_ids:
                    genotype = genotype_mapping[sid]
                    if not genotype in mutation_geno_counts[dna_name]:
                        mutation_geno_counts[dna_name][genotype]= 0
                    mutation_geno_counts[dna_name][genotype] += 1
    for mutation in mutation_geno_counts:
        for genotype in mutation_geno_counts[mutation]:
            perc = mutation_geno_counts[mutation][genotype] / genotype_sample_counts[genotype]
            if perc >= min_thresh and perc <= max_thresh:
                shared_mutations[genotype].append(mutation)


    return shared_mutations

def qa_genotype_kmers(genotype_mapping,kmer_geno_assoc,uid_map):
    report = {}

    # Get counts of each genotype
    genotype_sample_counts = {}
    for sample_id in genotype_mapping:
        genotype = genotype_mapping[sample_id]
        if not genotype in genotype_sample_counts:
            genotype_sample_counts[genotype] = 0
            report[genotype] = {
                'name': genotype,
                'num_members': 0,
                'num_shared_kmers': 0,
                'num_diagnostic_kmers': 0,
                'shared_kmers': [],
                'diagnostic_kmers': [],
                'qc_message': ''
            }
        genotype_sample_counts[genotype] += 1
        report[genotype]['num_members'] += 1

    for uid in kmer_geno_assoc:
        shared_genotypes = kmer_geno_assoc[uid]['shared']
        for genotype in shared_genotypes:
            n = uid_map[uid]
            report[genotype]['shared_kmers'].append(n)
        diagnostic_genotypes = kmer_geno_assoc[uid]['diagnostic']
        for genotype in diagnostic_genotypes:
            n = uid_map[uid]
            report[genotype]['diagnostic_kmers'].append(n)


    for genotype in report:
        report[genotype]['num_shared_kmers'] = len(report[genotype]['shared_kmers'])
        report[genotype]['num_diagnostic_kmers'] = len(report[genotype]['diagnostic_kmers'])
        report[genotype]['shared_kmers'].sort()
        report[genotype]['shared_kmers'] = ';'.join(report[genotype]['shared_kmers'] )
        report[genotype]['diagnostic_kmers'].sort()
        report[genotype]['diagnostic_kmers'] = ';'.join(report[genotype]['diagnostic_kmers'])

        qc_message = []
        if report[genotype]['num_shared_kmers'] == 0:
            qc_message.append( 'FAIL: No shared kmers' )
        if report[genotype]['diagnostic_kmers'] == 0:
            qc_message.append( 'WARNING: No diagnostic kmers' )
        report[genotype]['qc_message'] = ';'.join(qc_message)

    return report

def qa_genotype_mutations(genotype_mapping,mutation_geno_association):
    report = {}

    # Get counts of each genotype
    genotype_sample_counts = {}
    for sample_id in genotype_mapping:
        genotype = genotype_mapping[sample_id]
        if not genotype in genotype_sample_counts:
            genotype_sample_counts[genotype] = 0
            report[genotype] = {
                'name': genotype,
                'num_members': 0,
                'num_shared_mutations': 0,
                'num_diagnostic_mutations': 0,
                'shared_mutations': [],
                'diagnostic_mutations': [],
                'qc_message': ''
            }
        genotype_sample_counts[genotype] += 1
        report[genotype]['num_members'] += 1

    for mkey in mutation_geno_association:
        shared_genotypes = mutation_geno_association[mkey]['shared']
        for genotype in shared_genotypes:
            report[genotype]['shared_mutations'].append(mkey)
        diagnostic_genotypes = mutation_geno_association[mkey]['diagnostic']
        for genotype in diagnostic_genotypes:
            report[genotype]['diagnostic_mutations'].append(mkey)


    for genotype in report:
        report[genotype]['num_shared_mutations'] = len(report[genotype]['shared_mutations'])
        report[genotype]['num_diagnostic_mutations'] = len(report[genotype]['diagnostic_mutations'])
        report[genotype]['shared_mutations'].sort()
        report[genotype]['shared_mutations'] = ';'.join(report[genotype]['shared_mutations'] )
        report[genotype]['diagnostic_mutations'].sort()
        report[genotype]['diagnostic_mutations'] = ';'.join(report[genotype]['diagnostic_mutations'])

        qc_message = []
        if report[genotype]['num_shared_mutations'] == 0:
            qc_message.append( 'FAIL: No shared mutations' )
        if report[genotype]['diagnostic_mutations'] == 0:
            qc_message.append( 'WARNING: No diagnostic mutations' )
        report[genotype]['qc_message'] = ';'.join(qc_message)

    return report

def print_scheme(scheme,file,header):
    """
    Takes the kmer dictionary and writes it to a file
    :param scheme: dict of kmer info
    :param file: file path
    :return:
    """
    fh = open(file,'w')
    fh.write("\t".join(header) + "\n")
    for mutation_key in scheme:
        for state in scheme[mutation_key]:
            for data in scheme[mutation_key][state]:
                row = []
                for field in header:
                    if field in data:
                        row.append(data[field])
                    else:
                        row.append('')
                fh.write("{}\n".format("\t".join([str(x) for x in row])))
    fh.close()

def impute_ambiguous_bases(input_alignment,ref_id,min_frac=0.9):
    align_len = len(input_alignment[ref_id])
    consensus_bases = calc_consensus(input_alignment)
    num_seqs = len(input_alignment)
    for seq_id in input_alignment:
        if seq_id == ref_id:
            continue
        seq = list(input_alignment[seq_id])
        for i in range(0,align_len):
            base = seq[i]
            if base in ['A','T','C','G','-']:
                continue
            consensus_base = max(consensus_bases[i].items(), key=operator.itemgetter(1))[0]
            consensus_base_count = consensus_bases[i][consensus_base]
            if consensus_base in ['A','T','C','G'] and consensus_base_count/num_seqs > min_frac:
                seq[i] = consensus_base
        input_alignment[seq_id] = ''.join(seq)
    return input_alignment

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

def build_kmer_profiles(sequence_ids,scheme):
    '''
    :param scheme:
    :type scheme:
    :return:
    :rtype:
    '''
    profiles = {}
    mutations = list(scheme.keys())
    for sample_id in sequence_ids:
        profiles[sample_id] = {}
        for mutation_key in mutations:
            profiles[sample_id][mutation_key] = 0

    #only flip the alt status
    for mutation_key in scheme:
        for kmer_entry in scheme[mutation_key]['alt']:
            seq_ids = kmer_entry['seq_ids']
            for seq_id in seq_ids:
                profiles[seq_id][mutation_key] = 1

    return profiles

def multiplicity_check(scheme,input_alignment,genotype_mapping):
    kmers_to_uid = {}
    uid_to_kmer = {}
    uid_to_state = {}
    mut_to_uid = {}
    uid = 0
    for mutation in scheme:
        mut_to_uid[mutation] = []
        for state in ['ref','alt']:
            for i in range(0,len(scheme[mutation][state])):
                record = scheme[mutation][state][i]
                uid_to_kmer[uid] = record['unalign_kseq']
                uid_to_state[uid] = state
                if not record['unalign_kseq'] in kmers_to_uid:
                    kmers_to_uid[record['unalign_kseq']] = []
                kmers_to_uid[record['unalign_kseq']].append(uid)
                mut_to_uid[mutation].append(uid)
                record['key'] = uid
                uid+=1


    A = init_automaton_dict(uid_to_kmer)
    kmer_results = {}
    for sample_id in input_alignment:
        kmer_results[sample_id] = perform_kmerSearch_fasta(uid_to_kmer, kmers_to_uid, A, {sample_id:input_alignment[sample_id]},1)
        for kmer in kmers_to_uid:
            uids = kmers_to_uid[kmer]
            for uid in uids:
                if uid in kmer_results[sample_id]:
                    freq = kmer_results[sample_id][uid]
                    break
            for uid in uids:
                kmer_results[sample_id][uid] = freq





    info = {'warnings':{},'errors':{},'valid':{}}
    for sample_id in kmer_results:
        results = kmer_results[sample_id]
        for mutation in mut_to_uid:
            count_ref = 0
            count_alt = 0
            uids = mut_to_uid[mutation]

            for uid in uids:
                if results [uid] == 0:
                    continue
                state = uid_to_state[uid]
                if state == 'ref':
                    count_ref += results[uid]
                else:
                    count_alt += results[uid]

            if count_ref == count_alt and count_ref == 0:
                if not mutation in info['warnings']:
                    info['warnings'][mutation] = {}
                info['warnings'][mutation][sample_id] = 'missing'

            elif count_ref >= 1 and count_alt >= 1 :
                if not mutation in info['errors']:
                    info['errors'][mutation] = {}
                info['errors'][mutation][sample_id] = 'multi'

            else:
                if not mutation in info['valid']:
                    info['valid'][mutation] = {}
                info['valid'][mutation][sample_id] = 'valid'

    for mutation in scheme:
        for state in ['ref','alt']:
            for i in range(0,len(scheme[mutation][state])):
                if mutation in info['valid'] or mutation in info['errors']:
                    scheme[mutation][state][i]['is_kmer_found'] = True
                else:
                    scheme[mutation][state][i]['is_kmer_found'] = False
                    scheme[mutation][state][i]['is_valid'] = False
                if mutation in info['errors']:
                    scheme[mutation][state][i]['is_kmer_unique'] = False
                    scheme[mutation][state][i]['is_valid'] = False

    return scheme


def run():
    #get arguments
    cmd_args = parse_args()
    logger = init_console_logger(2)

    #input parameters
    input_msa = cmd_args.input_msa
    input_meta = cmd_args.input_meta
    prefix = cmd_args.prefix
    ref_id = cmd_args.ref_id
    ref_gbk = cmd_args.ref_gbk
    outdir = cmd_args.outdir
    min_len = cmd_args.min_len
    max_len = cmd_args.max_len
    max_ambig = cmd_args.max_ambig
    min_complexity = cmd_args.min_complexity
    max_missing = cmd_args.max_missing
    n_threads = cmd_args.n_threads
    impute_min_frac = cmd_args.iFrac
    jellyfish_path = os.path.join(outdir,"jellyfish")


    # initialize analysis directory
    if not os.path.isdir(outdir):
        logger.info("Creating analysis results directory {}".format(outdir))
        os.mkdir(outdir, 0o755)
    else:
        logger.info("Results directory {} already exits, will overwrite any results files here".format(outdir))


    #output files
    scheme_file = os.path.join(outdir,"{}-scheme.txt".format(prefix))
    genotypes_mut_file = os.path.join(outdir,"{}-mutation.associations.txt".format(prefix))
    genotypes_kmers_file = os.path.join(outdir, "{}-kmers.associations.txt".format(prefix))
    genotype_dendrogram = os.path.join(outdir, "{}-dendropgram.html".format(prefix))

    #Get the Gene features from the reference sequence
    ref_features = parse_reference_sequence(ref_gbk)

    #Read the input MSA and calculate the consensus sequence
    logger.info("Reading input alignment {}".format(input_msa))

    input_alignment = read_fasta(input_msa)
    if ref_id not in input_alignment:
        logger.error("Reference id specified '{}' is not found in alignment, make sure it is present and retry".format(ref_id))
        sys.exit()
    if cmd_args.impute:
        logger.info("Imputing ambiguous bases enabled with fraction of {}".format(impute_min_frac))
        input_alignment = impute_ambiguous_bases(input_alignment, ref_id, min_frac=impute_min_frac)
    else:
        logger.info("Imputing ambiguous bases not enabled")

    min_members = len(input_alignment) - len(input_alignment)*max_missing
    logger.info("Found {} sequences in msa".format(len(input_alignment)))
    consensus_bases = calc_consensus(input_alignment)
    consensus_seq = generate_consensus_seq(consensus_bases)

    #read the metadata associations
    logger.info("Reading genotype associations from {}".format(input_meta))
    metadata_df = read_tsv(input_meta)
    logger.info("Found {} lines in {}".format(len(metadata_df),input_meta))
    metadata_df['genotype'] = metadata_df['genotype'].astype(str)
    genotype_mapping = get_genotype_mapping(metadata_df)


    #Remove pseudo sequences generated by parsityper if applicable
    logger.info("Filtering samples from genotype map which are not also in the alignment")
    missing_samples = list(set(list(genotype_mapping.keys())) - set(list(input_alignment.keys())) )
    samples_to_mask = ['pseudo_A','pseudo_T','pseudo_C','pseudo_G','consensus']  + missing_samples
    logger.info("Masking {} samples".format(len(missing_samples)))
    for sample in samples_to_mask:
        if sample in genotype_mapping:
            logger.info("Removing sample {} from analysis due to no sequence in alignment".format(sample))
            del(genotype_mapping[sample])


    #Identify variable positions within the alignment
    logger.info("Scanning alignment for SNPs")
    snp_positions = find_snp_positions(consensus_seq)
    logger.info("Found {} variable sites".format(len(snp_positions)))
    sequence_deletions = {}
    logger.info("Scanning alignment for indels")
    for seq_id in input_alignment:
        sequence_deletions[seq_id] = find_internal_gaps(input_alignment[seq_id])
    variant_positions = []
    for pos in snp_positions:
        variant_positions.append([pos,pos])
    unique_indels = {}
    for sid in sequence_deletions:
        for i in sequence_deletions[sid]:
            unique_indels[i] = ''
    logger.info("Found {} potential indels".format(len(unique_indels)))
    for indel in unique_indels:
        (start,end) = indel.split(':')
        variant_positions.append([int(start),int(end)])

    logger.info("Creating scheme based on identified SNPs and indels")

    scheme = construct_scheme(variant_positions, ref_id, input_alignment, jellyfish_path, min_len, max_len,n_threads)

    logger.info("Performing QA on selected k-mers")
    scheme = qa_scheme(scheme, input_alignment, ref_id, genotype_mapping,min_len=min_len, max_len=max_len,min_complexity=min_complexity)

    logger.info("Adding gene annotations")
    scheme = add_gene_inference(scheme, ref_id, ref_features, input_alignment)

    #remove samples where genotype mapping info is missing
    valid_samples = set(list(genotype_mapping.keys()))
    unique_index = 0
    for mutation_key in scheme:
        for state in scheme[mutation_key]:
            for record in scheme[mutation_key][state]:
                record['seq_ids'] = list( set(record['seq_ids']) & valid_samples )
                record['key'] = unique_index
                unique_index+=1


    scheme = format_scheme_human_readable(scheme, ref_id, input_alignment)

    #get the kmer profile for each sample
    logger.info("Building genotype kmer profiles")
    kmer_profile = build_kmer_profiles(list(genotype_mapping.keys()),scheme)

    # create a plot of sample similarity for a multi-sample run
    if len(kmer_profile ) > 1:
        logger.info("Plotting Sample dendrogram")
        labels = []
        for sample in kmer_profile:
            if sample in genotype_mapping:
                genotype = genotype_mapping[sample]
            else:
                genotype = 'n/a'
            labels.append("{} | {}".format(sample,genotype))
        d = dendrogram_visualization()
        kmer_content_profile_df = pd.DataFrame.from_dict(kmer_profile,orient='index')
        d.build_dendrogram(labels, kmer_content_profile_df,
                                      genotype_dendrogram)
    #identify genotype shared kmers
    logger.info("Identifying shared kmers by genotype")
    shared_kmers = identify_shared_kmers(genotype_mapping,scheme,min_thresh=0.5,max_thresh=1)
    positive_kmers = identify_shared_kmers(genotype_mapping,scheme,min_thresh=0.95,max_thresh=1)
    partial_positive_kmers = identify_shared_kmers(genotype_mapping, scheme, min_thresh=0.01, max_thresh=0.94999)


    #identify shared mutations
    shared_mutations = identify_shared_mutations(genotype_mapping, scheme, min_thresh=0.5, max_thresh=1)
    positive_mutations = identify_shared_mutations(genotype_mapping, scheme, min_thresh=0.95,max_thresh=1)
    partial_positive_mutations = identify_shared_mutations(genotype_mapping, scheme, min_thresh=0.01, max_thresh=0.94999)

    #reorder datastructure to add diagnostic information to the scheme
    kmer_geno_assoc = {}

    for genotype in shared_kmers:
        for uid in shared_kmers[genotype]:
            if not uid in kmer_geno_assoc:
                kmer_geno_assoc[uid] = {'positive':[],'partial':[],'shared':[],'diagnostic':[]}
            kmer_geno_assoc[uid]['shared'].append(genotype)
        for uid in positive_kmers[genotype]:
            if not uid in kmer_geno_assoc:
                kmer_geno_assoc[uid] = {'positive':[],'partial':[],'shared':[],'diagnostic':[]}
            kmer_geno_assoc[uid]['positive'].append(genotype)
        for uid in partial_positive_kmers[genotype]:
            if not uid in kmer_geno_assoc:
                kmer_geno_assoc[uid] = {'positive':[],'partial':[],'shared':[],'diagnostic':[]}
            kmer_geno_assoc[uid]['partial'].append(genotype)

    uid_index = {}
    #Add positive and partial groupings to scheme
    for mutation_key in scheme:
        for state in scheme[mutation_key]:
            for row in scheme[mutation_key][state]:
                uid = row['key']
                uid_index[uid] = "{}.{}.{}".format(uid,row['dna_name'],row['state'])
                if uid in kmer_geno_assoc:
                    row['positive_genotypes'] = kmer_geno_assoc[uid]['positive']
                    row['partial_genotypes'] = kmer_geno_assoc[uid]['partial']

    mutation_geno_association = {}
    for genotype in shared_mutations:
        for mkey in shared_mutations[genotype]:
            if not mkey in mutation_geno_association:
                mutation_geno_association[mkey] = {'positive':[],'partial':[],'shared':[],'diagnostic':[]}
            mutation_geno_association[mkey]['shared'].append(genotype)
        for mkey in positive_mutations[genotype]:
            if not mkey in mutation_geno_association:
                mutation_geno_association[mkey] = {'positive':[],'partial':[],'shared':[],'diagnostic':[]}
            mutation_geno_association[mkey]['positive'].append(genotype)
        for mkey in partial_positive_mutations[genotype]:
            if not mkey in mutation_geno_association:
                mutation_geno_association[mkey] = {'positive':[],'partial':[],'shared':[],'diagnostic':[]}
            mutation_geno_association[mkey]['partial'].append(genotype)


    #identify genotype diagnostic kmers
    logger.info("Identifying diagnostic kmers by genotype")
    for uid in kmer_geno_assoc:
        shared = kmer_geno_assoc[uid]['shared']
        partial = kmer_geno_assoc[uid]['partial']
        if len(shared) == 1 or len(partial) == 1:
            genotype = shared[0]
            kmer_geno_assoc[uid]['diagnostic'].append(genotype)

    # identify genotype diagnostic mutations
    logger.info("Identifying diagnostic mutation events by genotype")
    for mkey in mutation_geno_association:
        shared = mutation_geno_association[mkey]['shared']
        partial = mutation_geno_association[mkey]['partial']
        if len(shared) == 1 or len(partial) == 1:
            genotype = shared[0]
            mutation_geno_association[mkey]['diagnostic'].append(genotype)

    #Analyse the genotypes for conflicts
    logger.info("Performing genotype kmer quality control analysis")
    genotype_report_kmer = qa_genotype_kmers(genotype_mapping,kmer_geno_assoc,uid_index)

    logger.info("Performing genotype mutation quality control analysis")
    genotype_report_mutation = qa_genotype_mutations(genotype_mapping,mutation_geno_association)


    #sort all multi-entry fields and collapse them to strings
    for mutation_key in scheme:
        for state in scheme[mutation_key]:
            for row in scheme[mutation_key][state]:
                row['mutation_key'] = mutation_key
                row['positive_genotypes'].sort()
                row['partial_genotypes'].sort()
                row['seq_ids'].sort()
                row['positive_genotypes'] = ','.join([str(x) for x in row['positive_genotypes']])
                row['partial_genotypes'] = ','.join([str(x) for x in row['partial_genotypes']])
                #deprecated negative TODO delete from templates
                row['negative_genotypes'] = ''
                row['seq_ids'] = ','.join([str(x) for x in row['seq_ids']])



    #write the final scheme to a file
    logger.info("Writting completed scheme")
    print_scheme(scheme, scheme_file,SCHEME_HEADER)


    # write the genotype file
    logger.info("Writting genotype report")
    pd.DataFrame.from_dict(genotype_report_mutation,orient='index').to_csv(genotypes_mut_file,sep='\t',index=False)
    pd.DataFrame.from_dict(genotype_report_kmer, orient='index').to_csv(genotypes_kmers_file, sep='\t', index=False)


run()


