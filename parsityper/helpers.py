from Bio import SeqIO
import pandas as pd
import logging, os, sys, re, collections, operator, math, shutil, datetime, copy, time
from itertools import product
from multiprocessing import Pool, Manager
from collections import Counter
from scipy.spatial.distance import cdist
from parsityper.constants import HTML_TEMPLATE_FILE, LOG_FORMAT, TYPING_SCHEMES, RUN_INFO, NEGATE_BASE_IUPAC, IUPAC_LOOK_UP, bases_dict
from parsityper.words import NOUNS, COLORS, DESCRIPTORS
from parsityper.kmerSearch import  init_automaton_dict,expand_degenerate_bases,find_in_fasta_dict,parallel_query_contigs
import random, hashlib
import numpy as np
from datetime import datetime
from Bio import GenBank
from Bio import SeqIO
from Bio.Seq import Seq
import glob
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from sklearn.metrics.pairwise import nan_euclidean_distances

NT_SUB = str.maketrans('acgtrymkswhbvdnxACGTRYMKSWHBVDNX',
                       'tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX')



def get_sequence_files(directory):
    '''
    Globs a directory and files sequence files
    :param directory: path
    :type directory: string
    :return: list of found files
    :rtype: list
    '''
    files = {'fasta':[],'fastq':[]}
    file_types = {
        '.fa':'fasta',
        '.fas':'fasta',
        'fasta':'fasta',
        '.fq': 'fastq',
        '.fastq': 'fastq',
    }
    for name in glob.glob(os.path.join(directory,"*")):
        for extension in file_types:
            if extension in name:
                file_type = file_types[extension]
                files[file_type].append(name)
                break
    return files




def expand_degenerate_bases(seq):
    """List all possible kmers for a scheme given a degenerate base

    Args:
         Scheme_kmers from SNV scheme fasta file

    Returns:
         List of all possible kmers given a degenerate base or not
    """

    return list(map("".join, product(*map(bases_dict.get, seq))))


def revcomp(s):
    """Reverse complement nucleotide sequence

    Args:
        s (str): nucleotide sequence

    Returns:
        str: reverse complement of `s` nucleotide sequence
    """
    return s.translate(NT_SUB)[::-1]

def get_expanded_kmer_number(scheme):
    """Checks that the number of kmers about to be created is not at too high a computation or time cost

    Args:
         scheme: kmer scheme dict

    """
    kmer_number = 0
    for kmer_id in scheme:
        seqs = [scheme[kmer_id]['positive'],scheme[kmer_id]['negative']]
        for seq in seqs:
            value = 1
            for char in seq:
                length_key = len(bases_dict[char])
                value *= length_key
            kmer_number += value
    return kmer_number


def generate_random_phrase():
    '''
    Generates a random phrase for naming kmer profiles
    '''
    phrase = []
    phrase.append(DESCRIPTORS[random.randrange(0, len(DESCRIPTORS) - 1)])
    phrase.append(COLORS[random.randrange(0, len(COLORS) - 1)])
    phrase.append(NOUNS[random.randrange(0,len(NOUNS)-1)])
    return(phrase)

def get_list_item(string,list_len):
    md5 = hashlib.md5()
    md5.update(string.encode())
    v = int(md5.hexdigest(), 16)
    r = v % list_len
    if r >= list_len -1:
        r = list_len -1

    return r

def generate_phase_md5(md5):
    '''
    Generates a random phrase for naming kmer profiles
    '''
    phrase = []
    list_len = len(DESCRIPTORS)
    i = get_list_item(md5, list_len)
    phrase.append(DESCRIPTORS[i])

    list_len = len(COLORS)
    i = get_list_item(md5, list_len)
    phrase.append(COLORS[i])

    list_len = len(NOUNS)
    i = get_list_item(md5, list_len)
    phrase.append(NOUNS[i])
    return(phrase)



def init_console_logger(lvl=2):
    '''
    Controlls the level of messaging provided to the user
    :param lvl: integer indicating level
    :return: logging object
    '''
    root = logging.getLogger()

    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]

    report_lvl = logging_levels[lvl]
    root.setLevel(report_lvl)  # set root logger level

    logging.basicConfig(format=LOG_FORMAT, level=report_lvl)

    return logging.getLogger(__name__)



def read_scheme_fasta(fasta_file):
    '''
    :param fasta_file: Fasta file
    :return: dictionary of sequences indexed by id
    '''
    reference = {}
    for seq_record in SeqIO.parse(fasta_file,format='fasta'):
        reference[str(seq_record.id)] = list(str(seq_record.seq))
    return reference

def report_template():
    '''
    Reads template file into string
    :return: HTML template string
    '''
    return open(HTML_TEMPLATE_FILE, 'r').read()

def read_tsv(tsv_file):
    '''
    Reads TSV file into a pandas data frame
    :param tsv_file: tsv file path
    :return: pd dataframe
    '''
    return pd.read_csv(tsv_file,sep="\t",header=0,low_memory=False)

def scheme_to_biohansel_fasta(df,fasta_file):
    '''
    Writes a compatible kmer df into a biohansel compatible scheme fasta
    :param df: kmer scheme data frame
    :param fasta_file: fasta file path
    :return:
    '''
    out_string = []
    for row in df.itertuples():
        id = row.key
        pos_kmer = row.positive
        neg_kmer = row.negative
        out_string.append(">{}-{}\n{}\n>negative{}-{}\n{}".format(id,id,pos_kmer,id,id,neg_kmer))
    fh = open(fasta_file,'w')
    fh.write("\n".join(out_string))
    fh.close

def filter_biohansel_kmer_df(df,min_freq,min_frac):
    df = df[df['freq'] > min_freq]
    df = df[df['pos_ratio'] > min_freq]
    return df

def calc_md5(string):
    seq = str(string).encode()
    md5 = hashlib.md5()
    md5.update(seq)
    return md5.hexdigest()

def validate_args(cmd_args,logger):
    '''
    Validates command line parameters
    :param cmd_args: args object
    :param logger: logging object
    :return: True on success, false if any parameter fails its check
    '''
    is_valid = True
    min_cov = cmd_args.min_cov
    if min_cov < 1:
        logger.error("Error you need to specify an integer >= 1 for min_cov")
        is_valid = False

    min_cov_frac = cmd_args.min_cov_frac
    if min_cov_frac <= 0 or min_cov_frac > 1:
        logger.error("Error you need to specify a float value >0 and <= 1")
        is_valid = False

    scheme= cmd_args.scheme
    if not scheme in TYPING_SCHEMES and not os.path.isfile(scheme):
        logger.error("Error specified scheme name or file does not exist".format(scheme))
        is_valid = False

    mode = cmd_args.mode
    if mode != 'single' and mode != 'batch':
        logger.error("Error specified operating mode is invalid, enter 'single' or 'batch'")
        is_valid = False

    type = cmd_args.type
    if type != 'single' and type != 'multi':
        logger.error("Error specified sample type is invalid, enter 'single' or 'multi' you entered {}".format(type))
        is_valid = False

    R1 = cmd_args.R1
    R2 = cmd_args.R2
    SE = cmd_args.se

    if R1 is not None and R2 is not None:
        if not os.path.isfile(R1):
            logger.error("Error {} is not found".format(R1))
            is_valid = False
        if not os.path.isfile(R2):
            logger.error("Error {} is not found".format(R2))
            is_valid = False

    if SE is not None:
        if not os.path.isfile(SE):
            logger.error("Error {} is not found".format(SE))
            is_valid = False

    if (R1 is not None or R2 is not None) and SE is not None:
        logger.error("Error you have specified both paired-end reads and single-end, please specify only one")
        is_valid = False

    data_dir = cmd_args.data_dir
    if data_dir is not None:
        if not os.path.isdir(data_dir):
            logger.error("Error {} directory is not found".format(data_dir))
            is_valid = False

    if ((R1 is not None or R2 is not None) or SE is not None) and data_dir is not None:
        logger.error("Error you have specified readsets and a data directory, please specify only one")
        is_valid = False

    if (R1 is None or R2 is None) and SE is None and data_dir is None:
        logger.error("Error you need to specify either read sets or a data directory to process")
        is_valid = False

    return is_valid


def init_kmer_targets(scheme_df):
    '''

    :param scheme_df: K-mer typing scheme pandas dataframe
    :return: dict of kmer targets
    '''
    targets = {}
    fields = list(set(scheme_df.columns.tolist()) )
    for row in scheme_df.itertuples():
        key = str(row.key)
        if not key in targets:
            targets[key] = {}
        for field in fields:
            targets[key][field] = getattr(row, field)

    return targets


def generate_biohansel_kmer_names(value_list):
    '''

    :param value_list: integer list
    :return: list of kmer names in biohansel format
    '''
    names = []
    for value in value_list:
        names.append("{}-{}".format(value,value))
        names.append("negative{}-{}".format(value, value))
    return names

def get_scheme_template(elements,value):
    template = dict()
    for e in elements:
        if isinstance(value,dict) or isinstance(value,list) :
            value = copy.deepcopy(value)
        template[e] = value
    return template

def calc_profile_distances(sample_profile,genotype_profiles):
    '''
    Computes the distance between the sample profile and all of the precomputed genotypes
    :param sample_profile: pandas df of kmer profile, 0,0.5,1 or nan as values for sample
    :param genotype_profiles: pandas df of kmer profile, 0,0.5,1 or nan as values for reference to compare sample
    :return: list of distances
    '''
    return cdist(sample_profile, genotype_profiles, 'euclidean')[0]

def identify_candidate_genotypes_by_dist(sample_profile,genotype_profiles,max_dist=0.5):
    '''
    Computes the distance between the sample profile and all of the precomputed genotypes and filteres out any whihc
    are too far from the distance threshold
    :param sample_profile: pandas df of kmer profile, 0,0.5,1 or nan as values for sample
    :param genotype_profiles: pandas df of kmer profile, 0,0.5,1 or nan as values for reference to compare sample
    :return: list of distances
    '''
    dists = calc_profile_distances(sample_profile,genotype_profiles)
    samples = genotype_profiles.index.tolist()
    candidates = {}
    for i in range(0,len(dists)):
        if dists[i] > max_dist:
            continue
        candidates[samples[i]] = dists[i]
    return candidates

def profile_pairwise_distmatrix(profile_st):
    '''
    Computes pairwise jaccard distances between sample profiles
    '''
    samples = list(profile_st.keys())
    num_samples = len(samples)
    matrix = np.zeros((num_samples, num_samples))
    for i in range(0,len(samples)):
        for k in range(i,len(samples)):
            if i == k:
                jaccard = 0
            elif len(profile_st[samples[i]]) > 0 or len(profile_st[samples[k]]) >0:
                jaccard = 1 - (len(list(set(profile_st[samples[i]])  &  set(profile_st[samples[k]]))) / \
                          len(list(set(profile_st[samples[i]]) | set(profile_st[samples[k]]))))
            else:
                jaccard = 1

            matrix[i,k] = jaccard
            matrix[k, i] = jaccard

    return matrix

def generate_run_report(sample_report,
                        report_sample_composition_summary,
                        sample_kmer_data, scheme_kmer_target_info,
                        sample_complexity,
                        max_mixed_sites, min_cov_frac,  min_cov,primers):
    RUN_INFO = {
        'Analysis Date': '',
        'Primer set': '',
        'Min coverage': 0,
        'Min coverage fraction': 0,
        'Input data': '',
        'Number of Samples': 0,
        'Average number of kmers per sample': 0,
        'Stdev number of kmers per sample': 0,
        'Average kmer coverage per sample': 0,
        'Stdev kmer coverage per sample': 0,
        'Average mixed kmers per sample': 0,
        'Stdev mixed kmers per sample': 0,
        'Number of QC Pass samples': 0,
        'Number of QC Warning samples': 0,
        'Number of QC Fail samples': 0,
        'Number of unique genotypes': 0,
        'Number of distinct profiles': 0,
    }
    run_info = RUN_INFO
    run_info['Analysis Date'] = datetime.today().strftime('%Y-%m-%d')
    run_info['Min coverage'] = min_cov
    run_info['Min coverage fraction'] = min_cov_frac
    run_info['Primer set'] = primers
    run_info['Number of Samples'] = len(sample_report)

    return run_info


def parse_reference_sequence(gbk_file):
    """
    :param gbk_file: Reference genbank format file with sequence annotations
    :return: dict of all of the reference features
    """
    with open(gbk_file) as handle:
        for record in GenBank.parse(handle):
            gb_accession = record.accession[0]
            gb_accession_version = record.version.split('.')[1]
            genome_seq = repr(record.sequence).replace("\'",'')
            sequences = {}
            sequences[gb_accession] = {
                'accession':gb_accession,
                'version': gb_accession_version,
                'features': {'source': genome_seq}
            }
            features = record.features
            for feat in features:
                if feat.key == 'CDS' or feat.key == '5\'UTR' or feat.key == '3\'UTR':
                    if not feat.key in sequences[gb_accession]['features']:
                        sequences[gb_accession]['features'][feat.key] = []
                    qualifier = feat.qualifiers
                    positions = []
                    gene_name = ''
                    aa = ''
                    for name in qualifier:
                        if name.key == '/gene=':
                            gene_name = name.value.replace("\"", '').strip()
                        if name.key == '/translation=':
                            aa = name.value.replace("\"", '').strip()

                    locations = feat.location.strip().replace("join(", '').replace(')', '').split(',')
                    seq = []

                    for location in locations:
                        location = location.replace('<','').replace('>','')
                        if not 'complement' in location:
                            location = location.split('.')
                            start = int(location[0]) - 1
                            end = int(location[2])
                            seq.append(genome_seq[start:end].replace("\'", ''))
                            positions.append([start, end])
                        else:
                            location = location.replace('complement(','').replace(')','').split('.')
                            start = int(location[0]) - 1
                            end = int(location[2])
                            seq.append(revcomp(genome_seq[start:end].replace("\'", '')))
                            positions.append([start, end])

                    seq = ''.join(seq)
                    sequences[gb_accession]['features'][feat.key].append(
                        {'gene_name': gene_name, 'dna_seq': seq, 'aa_seq': aa, 'positions': positions})

    return sequences



def get_aa_delta(start, end, variant, ref_info,ref_name,trans_table=1):
    #print("start {} end {}".format(start,end))
    vlen = end - start
    is_cds = False
    gene = ''
    gene_start = -1
    gene_end = -1
    aa_start = -1
    aa_end = -1
    ref_seq = ''
    alt_seq = ''
    ref_target = []
    alt_target = []
    is_frame_shift = False

    cds_start = start
    cds_end = cds_start + vlen

    is_silent = True


    count_gaps = variant.count('-')
    if count_gaps % 3 > 0:
        is_frame_shift = True

    found = False
    for feat in ref_info[ref_name]['features']['CDS']:
        positions = feat['positions']
        dna_seq = feat['dna_seq']
        aa_seq = str(Seq(feat['dna_seq']).translate(trans_table))
        gene_start = -1
        gene_end = -1
        spacer = 0
        for s, e in positions:
            if gene_start == -1:
                gene_start = s
            if gene_end < e:
                gene_end = e
            if cds_start  >= s and cds_end  <= e:

                gene = feat['gene_name']
                cds_start -= gene_start
                r = int(cds_start % 3)
                cds_start -= r
                spacer = r
                vlen += r
                r = vlen % 3
                num_codons = int((vlen - r) / 3)
                cds_end = cds_start + (num_codons * 3) + 1
                length = cds_end - cds_start
                r = int(length % 3)
                cds_end += r +1
                num_codons = int((cds_end - cds_start ) / 3)
                is_cds = True
                found = True
                ref_target = list(dna_seq[cds_start:cds_end])
                alt_target = list(dna_seq[cds_start:cds_end])
                aa_start = int(cds_start / 3)
                aa_end = aa_start + num_codons - 1
                break
        if found:
            break

    if len(ref_target) > 0:
       #mutate sequence to match variant
        #print(alt_target)
        for i in range(0,len(variant)):
            #if i + spacer >= len(alt_target):
                #print("{} out of range {} {} {}".format(i,len(variant),variant,alt_target))
            alt_target[i+spacer] = variant[i]

        #print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(start, end, cds_start, cds_end, variant, vlen, ref_target,alt_target))

        ref_seq = "{}".format(Seq(''.join(ref_target).replace('-','N')).translate(table=trans_table))
        alt_seq = "{}".format(Seq(''.join(alt_target).replace('-','N')).translate(table=trans_table))

        if ref_seq != alt_seq:
            is_silent = False
    gene_start += 1
    gene_end += 1
    aa_start += 1
    aa_end += 1
    if not is_cds:
        ref_seq = ''
        alt_seq = ''
    if len(gene) == 0:
        gene_start = -1
        gene_end = -1
        aa_start = -1
        aa_end = -1
        cds_start = -1
        cds_end = -1
    return {
        'gene': gene,
        'gene_start':gene_start,
        'gene_end':gene_end,
        'cds_start':cds_start,
        'cds_end':cds_end,
        'aa_start':aa_start,
        'aa_end':aa_end,
        'ref_state':ref_seq,
        'alt_state':alt_seq,
        'is_silent':is_silent,
        'is_cds': is_cds,
        'is_frame_shift':is_frame_shift
    }


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


def calc_consensus(input_alignment):
    """
    Determines the character counts for each position in an MSA
    :param input_alignment: dict of sequences
    :return: dict of base counts by position
    """
    seq_id = next(iter(input_alignment))
    seq = input_alignment[seq_id]
    seq_len = len(seq)

    # get consensus
    consensus = []
    for i in range(0, seq_len):
        consensus.append({'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0, '-': 0})

    for seq_id in input_alignment:
        seq = input_alignment[seq_id]
        #if seq_len != len(seq):
        #    print(seq_id)
        for i in range(0,seq_len):
            base = seq[i].upper()
            if base in consensus[i]:
                consensus[i][base]+=1
            else:
                consensus[i]['N'] += 1

    return consensus

def generate_consensus_seq(consensus):
    """
    Using the base frequency from the multiple sequence alignment, it generates a pseudo consensus sequence
    with variable positions masked with IUPAC characters
    :param consensus: dict of base counts by position
    :return: str DNA consensus sequence
    """
    cons_seq = []
    bases = ['A','T','C','G']
    iupac = IUPAC_LOOK_UP
    variable_sites = []

    for i in range(0,len(consensus)):
        lookup_key = []
        sum = 0
        for base in consensus[i]:
            sum+= consensus[i][base]
        for base in consensus[i]:
            if base in bases and consensus[i][base] > 0:
                lookup_key.append(base)
        if consensus[i]['-'] / sum > 0.9:
            lookup_key = '-'
        else:
            lookup_key = sorted(lookup_key)
            lookup_key = ''.join(lookup_key)
        if lookup_key in iupac:
            base = iupac[lookup_key]
        else:
            base = max(consensus[i].items(), key=operator.itemgetter(1))[0]
        if base not in ['A','T','C','G','-']:
            variable_sites.append(i)
        cons_seq.append(base)
    return ''.join(cons_seq)

def read_fasta(fasta_file):
    """
    Reads fasta file into a dict
    :param fasta_file: fasta formatted file
    :return: dict of sequences
    """
    reference = {}
    for seq_record in SeqIO.parse(fasta_file,format='fasta'):
        reference[str(seq_record.id)] = str(seq_record.seq).upper()
    return reference

def find_gaps(seq):
    """
    Accepts a string and returns the positions of all of the gaps in the sequence
    :param seq: str
    :return: list of [start,end] of all of the gaps
    """
    match = re.finditer(r"-+", seq)
    positions = []
    for m in match:
        positions.append([m.start(),m.end()])
    return positions

def find_internal_gaps(seq):
    """
    Accepts a string and returns the positions of all of the gaps in the sequence which are flanked by nt bases
    :param seq: str
    :return: list of [start,end] of all of the internal gaps
    """
    gaps = find_gaps(seq)
    seq_len =  len(seq) -1
    internal_gaps = []
    iupac = IUPAC_LOOK_UP
    for gap in gaps:
        start = gap[0]
        end = gap[1]
        if start == 0 or end >= seq_len:
            continue

        if seq[start-1] not in iupac or seq[end+1] not in iupac:
            continue
        internal_gaps.append("{}:{}".format(start,end))
    return internal_gaps

def find_snp_positions(consensus):
    """
    Acccepts IUPAC masked pseudo consensus sequence and identifies the positions which are not A,T,C,G
    :param consensus: str sequence
    :return: list of variable positions
    """
    positions = []
    for i in range(0,len(consensus)):
        base = consensus[i]
        if base not in ['A','T','C','G','N','-']:
            positions.append(i)
    return positions

def get_kmers(start,end,input_alignment):
    """
    Accepts as start and end position within a MSA and returns a dict of all of the kmers
    :param start: int
    :param end: int
    :param input_alignment: dict of sequences
    :return: dict of sequence kmers corresponding to the positions
    """
    kmers = {}
    seq_iter = input_alignment.keys()
    for seq_id in seq_iter :
        kmers[seq_id] = input_alignment[seq_id][start:end]
    return kmers


def find_initial_start(pos,reference_sequence,min_length):
    """
    Using an initial position, it finds the initial starting position which satisfies the minimum length
    :param pos: int
    :param reference_sequence: str
    :param min_length: int
    :return: int
    """
    ref_len=len(reference_sequence)
    nt_count = 0
    start = pos - 1
    for i in range(0,ref_len):
        base = reference_sequence[start].upper()
        if base != '-':
            nt_count += 1
        if start < 0:
            start = 0
            break
        if nt_count >= min_length:
            break
        start -= 1
    return start

def find_initial_end(pos,reference_sequence,min_length):
    """
    Using an initial position, it finds the initial ending position which satisfies the minimum length
    :param pos: int
    :param reference_sequence: str
    :param min_length: int
    :return: int
    """
    ref_len=len(reference_sequence)
    nt_count = 0
    end = pos - 1

    for i in range(0,ref_len):
        if end < 0 or end >= ref_len:
            end = ref_len
            break
        base = reference_sequence[end].upper()
        if base in ['A', 'T', 'C', 'G']:
            nt_count += 1

        if nt_count >= min_length:
            break
        end += 1
    return end


def count_kmers(seq, K=2):
    mers = {}
    """Count kmers in sequence"""
    for i in range(0,len(seq)):
        mer = list(seq[i:i+K])
        mer.sort()
        mer = ''.join(mer)
        if len(mer) != K:
            continue
        if not mer in mers:
            mers[mer] = 0
        mers[mer]+=1
    return mers


def optimize_kmer(pos,aln_seqs,reference_sequence,min_length,max_length,min_members,max_ambig=5,min_complexity=0.2,n_threads=1):
    """
    Accepts a position and a sequence and determines the kmer stretch which maximizes length, complexity and minimizes
    ambiguous characters
    :param pos: int position
    :param reference_sequence: str reference sequence
    :param min_length: int minimum length of kmer
    :param max_length: int maximum length of kmer
    :param max_ambig:  int maximum number of iupac characters
    :param min_complexity: float maximum percentage composition of one 2-mer
    :return:
    """

    prev_score = 0
    opt_kmer = [-1,-1]
    rlen = len(reference_sequence)
    if pos > 0:
        start = find_initial_start(pos, reference_sequence, min_length)
        # set fall back kmer
        istart = start
        iend = pos + 1
    else:
        start = 0
        istart = 0
        iend = min_length+1

    if iend > rlen:
        iend = rlen - 1
    if istart < 0:
        istart = 0
    kmers = {}
    positions = {}

    for length_target in range(min_length,max_length):
        for k in range(start ,pos):
            s = pos - k
            if s > length_target :
                continue
            rel_start = k
            nt_count = 0
            rel_end = k
            base_count = 0
            while nt_count < length_target:
                if base_count >= length_target or rel_end >= rlen -1:
                    break
                base = reference_sequence[rel_end]
                if base in ['A', 'T', 'C', 'G']:
                    nt_count += 1
                if base != '-':
                    base_count+=1
                rel_end += 1

            if start < 0 or start >= rlen or rel_end >= rlen or rel_end < pos:
                continue
            kmer = reference_sequence[rel_start:rel_end].replace('-','')
            klen = len(kmer)


            if klen > max_length :
                continue

            #count ambiguous characters
            bases = ['A','T','C','G']
            nt_count = 0
            mers = count_kmers(kmer, K=1)
            for b in bases:
                if b in mers:
                    nt_count+=mers[b]
            count_ambig = klen - nt_count

            #determine the complexity of the sequence and remove kmers composed heavily of the same 2-mer
            mers = count_kmers(kmer, K=2)
            num_mers = sum(mers.values())
            mer_perc = []

            for m in mers:
                mer_perc.append(mers[m]/num_mers )

            if len(mer_perc) > 0:
                minimum = min(mer_perc)
            else:
                minimum = 1
            score = (1 - ((nt_count+count_ambig)/max_length)) + (1 - minimum) + (1-count_ambig/max_length)
            if count_ambig <= max_ambig:
                kmers[kmer] = score
                positions[kmer] = [rel_start,rel_end]


    kmers = sorted(kmers.items(), key=lambda x: x[1], reverse=True)

    for kmer,score in kmers:

        count_ambig = len(kmer) - (kmer.count('A') + kmer.count('T') + kmer.count('C') + kmer.count('G'))
        if count_ambig > max_ambig or len(kmer) < min_length:
            continue
        A = init_automaton_dict({kmer:kmer})
        df = find_in_fasta_dict(A, aln_seqs)
        if is_kmer_valid(df, min_members):
            opt_kmer = opt_kmer = [positions[kmer][0], positions[kmer][1]]
            break

    if opt_kmer[0] == -1:
        opt_kmer = [istart, iend]
    kmer = reference_sequence[opt_kmer[0]:opt_kmer[1]].replace('-', '')

    if len(kmer) < min_length or len(kmer) > max_length:
        opt_kmer = [istart, iend]
    return opt_kmer

def add_snp_kmer_to_scheme(pos,ref_len,input_alignment,consensus_bases,consensus_seq,reference_info,ref_name,min_len,max_len,max_ambig,min_members,min_complexity=0.6,n_threads=1):
    scheme = {}
    from os import getpid
    stime = time.time()
    anything_but_bases = NEGATE_BASE_IUPAC
    ref_non_gap_lookup = generate_non_gap_position_lookup(input_alignment[ref_name])
    scheme = {}
    # remove any cases of variable site where it is a IUPAC character as the delta
    bases = consensus_bases[pos]
    snps = []
    for base in bases:
        if base in ['A', 'T', 'C', 'G']:
            if bases[base] > 0:
                snps.append(base)
    count_states = len(snps)
    if count_states == 1:
        return {}

    (start, end) = optimize_kmer(pos=pos,
                                 aln_seqs=input_alignment,
                                 reference_sequence=consensus_seq,
                                 min_length=min_len,
                                 max_length=max_len,
                                 min_members=min_members,
                                 max_ambig=max_ambig,
                                 min_complexity=min_complexity,
                                 n_threads=n_threads)

    if start < 0:
        start = 0

    rel_start = start
    rel_end = end + 1

    kmer = consensus_seq[start:end].replace('-', '')

    for i in range(0, len(snps)):
        base = snps[i]

        rel_pos = pos - rel_start
        pos_kmer = list(consensus_seq[rel_start:rel_end])
        pos_kmer[rel_pos] = base
        neg_kmer = list(consensus_seq[rel_start:rel_end])
        neg_kmer[rel_pos] = anything_but_bases[base]
        ref_base = input_alignment[ref_name][pos]
        if ref_base == base:
            continue

        non_gap_pos = ref_non_gap_lookup[pos]
        while non_gap_pos == -1:
            pos -= 1
            non_gap_pos = ref_non_gap_lookup[pos]

        kmer_name = "{}{}{}".format(base, non_gap_pos+1, ref_base)
        aa_info = get_aa_delta(non_gap_pos, non_gap_pos, base, reference_info, ref_name, trans_table=1)

        if len(aa_info['gene']) != 0 :
            kmer_name_aa = "{}{}{}".format(aa_info['ref_state'], aa_info['aa_start'], aa_info['alt_state'])
        else:
            kmer_name_aa = 'N/A'

        scheme[kmer_name] = {
            'dna_name': kmer_name,
            'aa_name': kmer_name_aa,
            'type': 'snp',
            'gene': aa_info['gene'],
            'gene_start': aa_info['gene_start'],
            'gene_end': aa_info['gene_end'],
            'cds_start': aa_info['cds_start']+1,
            'cds_end': aa_info['cds_end'],
            'is_silent': aa_info['is_silent'],
            'is_frame_shift': False,
            'ref_aa': aa_info['ref_state'],
            'alt_aa': aa_info['alt_state'],
            'variant_start': non_gap_pos+1,
            'variant_end': non_gap_pos+1,
            'kmer_start': rel_start + 1,
            'kmer_end': rel_end + 1,
            'variant_pos': base,
            'variant_neg': anything_but_bases[base],
            'positive': ''.join(pos_kmer),
            'negative': ''.join(neg_kmer),
            'positive_seqs': [],
            'partial_positive_seqs': [],
        }
        if aa_info['cds_start'] == -1:
            scheme[kmer_name]['cds_start'] = -1
        seq_bases = get_kmers(pos, pos + 1, input_alignment)
        for seq_id in seq_bases:
            if seq_bases[seq_id] == base:
                scheme[kmer_name]['positive_seqs'].append(seq_id)
        if len(scheme[kmer_name]['positive_seqs']) == 0:
            del (scheme[kmer_name])

    return scheme


def find_snp_kmers(input_alignment,snp_positions,consensus_bases,consensus_seq,reference_info,ref_name,min_len,max_len,max_ambig,min_members,min_complexity=0.6,n_threads=1):
    scheme = {}
    ref_len = len(input_alignment[ref_name])
    res = []
    pool = Pool(processes=n_threads)

    for pos in snp_positions:
        res.append(pool.apply_async(add_snp_kmer_to_scheme, args=(pos, ref_len, input_alignment, consensus_bases,
                                           consensus_seq, reference_info, ref_name,
                               min_len, max_len, max_ambig, min_members, 0.6, 1)))
    pool.close()
    pool.join()
    for x in res:
        if x._success:
            scheme.update(x.get())

    return scheme

def find_snp_kmers_bck(input_alignment,snp_positions,consensus_bases,consensus_seq,reference_info,ref_name,min_len,max_len,max_ambig,min_members,min_complexity=0.6,n_threads=1):
    scheme = {}
    anything_but_bases = NEGATE_BASE_IUPAC
    ref_len = len(input_alignment[ref_name])
    ref_non_gap_lookup = generate_non_gap_position_lookup(input_alignment[ref_name])
    used_kmer_positions = []

    # Add snps into the kmer scheme
    for pos in snp_positions:
        stime = time.time()
        scheme.update(add_snp_kmer_to_scheme(pos, ref_len, input_alignment, consensus_bases, consensus_seq, reference_info, ref_name,
                               min_len, max_len, max_ambig, min_members, min_complexity=0.6, n_threads=1))

        (start,end) = optimize_kmer(pos=pos,
                                    aln_seqs=input_alignment,
                                    reference_sequence=consensus_seq,
                                    min_length=min_len,
                                    max_length=max_len,
                                    min_members=min_members,
                                    max_ambig=max_ambig,
                                    min_complexity=min_complexity,
                                    n_threads=n_threads)
        kmer = consensus_seq[start:end].replace('-','')
        #print("{}\t{}\t{}".format(start,end,kmer))

        #print(time.time() - stime)
        if start < 0:
            start = 0

        rel_start = start
        rel_end = end + 1
        bases = consensus_bases[pos]

        #remove any cases of variable site where it is a IUPAC character as the delta
        snps = []
        for base in bases:
            if base in ['A', 'T', 'C', 'G']:
                if bases[base] > 0:
                    snps.append(base)
        count_states = len(snps)

        if count_states == 1:
            continue

        for i in range(0, len(snps)):
            base = snps[i]
            if i > 0:
                is_used = False
                while is_used:
                    is_used = False
                    for s, e in used_kmer_positions:
                        if s == rel_start and e == rel_end:
                            is_used = True
                            if e == rel_end:
                                rel_end += 1
                            if s == rel_start:
                                rel_start += 1
                            break
                    if rel_start > pos or rel_end > ref_len:
                        break

                if rel_start > pos or rel_end > ref_len:
                    rel_start = start
                    rel_end = end + 1


            rel_pos = pos - rel_start
            pos_kmer = list(consensus_seq[rel_start:rel_end])
            pos_kmer[rel_pos] = base
            neg_kmer = list(consensus_seq[rel_start:rel_end])
            neg_kmer[rel_pos] = anything_but_bases[base]
            ref_base = input_alignment[ref_name][pos]
            if ref_base == base:
                continue

            non_gap_pos = ref_non_gap_lookup[pos]
            while non_gap_pos == -1:

                pos -= 1
                non_gap_pos = ref_non_gap_lookup[pos]

            kmer_name = "{}{}{}".format(base, non_gap_pos+1, ref_base)
            aa_info = get_aa_delta(non_gap_pos, non_gap_pos, base, reference_info,ref_name,trans_table=1)

            if len(aa_info['gene']) != 0:
                kmer_name_aa = "{}{}{}".format(aa_info['ref_state'], aa_info['aa_start'], aa_info['alt_state'])
            else:
                kmer_name_aa = 'N/A'

            scheme[kmer_name] = {
                'dna_name': kmer_name,
                'aa_name': kmer_name_aa,
                'type': 'snp',
                'gene':aa_info['gene'],
                'gene_start':aa_info['gene_start'],
                'gene_end': aa_info['gene_end'],
                'cds_start': aa_info['cds_start']+1,
                'cds_end': aa_info['cds_end'],
                'is_silent':aa_info['is_silent'],
                'is_frame_shift': False,
                'ref_aa':aa_info['ref_state'],
                'alt_aa': aa_info['alt_state'],
                'variant_start': non_gap_pos+1,
                'variant_end': non_gap_pos+1,
                'kmer_start': rel_start + 1,
                'kmer_end': rel_end + 1,
                'variant_pos': base,
                'variant_neg': anything_but_bases[base],
                'positive': ''.join(pos_kmer),
                'negative': ''.join(neg_kmer),
                'positive_seqs': [],
                'partial_positive_seqs': [],

            }

            seq_bases = get_kmers(pos, pos + 1, input_alignment)
            for seq_id in seq_bases:
                if seq_bases[seq_id] == base:
                    scheme[kmer_name]['positive_seqs'].append(seq_id)
            if len(scheme[kmer_name]['positive_seqs']) == 0:
                del (scheme[kmer_name])

    return scheme


def add_indel_kmer_to_scheme():
    scheme = {}

    for indel in indels[seq_id]:
        (vstart,vend) = indel.split(':')
        vstart = int(vstart)
        vend = int(vend)
        vlen = (vend - vstart)

       #determine variant type
        ref_state = input_alignment[ref_name][vstart:vend +1].replace('-','')
        alt_state = input_alignment[seq_id][vstart:vend +1].replace('-','')

        if ref_state == alt_state:
            continue

        if len(ref_state) > len(alt_state):
            type = 'del'
        else:
            type = 'ins'

        if type == 'ins':
            (start, end) = optimize_kmer(pos=vend,
                                         aln_seqs=input_alignment,
                                         reference_sequence=consensus_seq,
                                         min_length=min_len,
                                         max_length=max_len,
                                         min_members=min_members,
                                         max_ambig=max_ambig,
                                         min_complexity=min_complexity,
                                         n_threads=n_threads)


        else:
            #(start, end) = optimize_kmer(vend, input_alignment[seq_id], min_len, max_len, max_ambig, min_complexity)
            (start, end) = optimize_kmer(pos=vend,
                                         aln_seqs=input_alignment,
                                         reference_sequence=input_alignment[seq_id],
                                         min_length=min_len,
                                         max_length=max_len,
                                         min_members=min_members,
                                         max_ambig=max_ambig,
                                         min_complexity=min_complexity)

        if (start == -1 or end == -1) :
            continue

        if vstart < start:
            if type == 'ins':
                #params = optimize_kmer(vstart, consensus_seq, min_len, max_len, max_ambig, min_complexity=0.5)

                params = optimize_kmer(pos=vstart,
                                             aln_seqs=input_alignment,
                                             reference_sequence=consensus_seq,
                                             min_length=min_len,
                                             max_length=max_len,
                                             min_members=min_members,
                                             max_ambig=max_ambig,
                                             min_complexity=min_complexity)

            else:
                #params = optimize_kmer(vstart, input_alignment[seq_id], min_len, max_len, max_ambig,
                #                             min_complexity=0.5)

                params = optimize_kmer(pos=vstart,
                                             aln_seqs=input_alignment,
                                             reference_sequence=input_alignment[seq_id],
                                             min_length=min_len,
                                             max_length=max_len,
                                             min_members=min_members,
                                             max_ambig=max_ambig,
                                             min_complexity=min_complexity)
            start = params[0]

        rel_start = (vstart - start)
        rel_end = rel_start + vlen

        if rel_end > end:
            rel_end = params[1]

        if (start == -1 or end == -1) :
            continue

        neg_kmer = list(consensus_seq[start:end+1])
        pos_kmer = list(consensus_seq[start:end+1])

        variant_pos = []
        variant_neg = []

        if type == 'ins':
            for i in range(rel_start,rel_end):

                neg_kmer[i] = '-'
                variant_neg.append('-')
                variant_pos.append(pos_kmer[i])

        if type == 'del':
            for i in range(rel_start,rel_end):

                pos_kmer[i] = '-'
                variant_neg.append(neg_kmer[i])
                variant_pos.append('-')

        variant_pos = ''.join(variant_pos)
        variant_neg = ''.join(variant_neg)

        #Handle long indels which exceed the max length of the kmer
        pos_kmer = ''.join(pos_kmer)
        neg_kmer = ''.join(neg_kmer)
        neg_kmer = neg_kmer.replace('-','')
        pos_kmer = pos_kmer.replace('-', '')
        if len(neg_kmer) > max_len:
            diff = False
            for i in range(0,len(neg_kmer)):
                if i >= len(pos_kmer):
                    break
                if neg_kmer[i] != pos_kmer[i]:
                    diff = True
                    end = start + i
                if diff and i >= max_len:
                    break
            if (end - start) < min_len:
                end = start + min_len
                i = min_len

        neg_kmer = neg_kmer[0:i+1]
        pos_kmer = pos_kmer[0:i+1]

        if len(neg_kmer) > min_len and len(pos_kmer) > min_len:
            #Trim any degenerate sites at the end of the sequence as long as that isn't a difference
            if neg_kmer[0] == pos_kmer[0]:
                if neg_kmer[0] not in ['A','T','C','G']:
                    start-=1
                    neg_kmer = neg_kmer[1:len(neg_kmer)]
                    pos_kmer = pos_kmer[1:len(neg_kmer)]
            if len(neg_kmer) == len(pos_kmer):
                if neg_kmer[len(neg_kmer)-1] == pos_kmer[len(pos_kmer)-1]:
                    if neg_kmer[len(neg_kmer)-1] not in ['A','T','C','G']:
                        end -=1
                        neg_kmer = neg_kmer[0:len(neg_kmer)-1]
                        pos_kmer = pos_kmer[0:len(pos_kmer)-1]

        neg_kmer = list(neg_kmer)
        pos_kmer = list(pos_kmer)
        if variant_pos == variant_neg:
            continue

        non_gap_start = ref_non_gap_lookup[vstart]
        pos = vstart
        while non_gap_start == -1:
            pos -= 1
            non_gap_start = ref_non_gap_lookup[pos]

        non_gap_end  = non_gap_start + vlen


        kmer_name = "{}{}_{}".format(type,non_gap_start, non_gap_end)
        aa_info = get_aa_delta(non_gap_start, non_gap_end, variant_pos, reference_info,ref_name,trans_table=1)

        if aa_info['aa_start'] != -1 and aa_info['aa_end'] !=-1:
            kmer_name_aa = "{}{}{}".format(aa_info['ref_state'], aa_info['aa_start'], aa_info['alt_state'])
        else:
            kmer_name_aa = 'N/A'
        if kmer_name not in scheme:
            scheme[kmer_name] = {
                'dna_name': kmer_name,
                'aa_name': kmer_name_aa,
                'type': type,
                'gene':aa_info['gene'],
                'gene_start':aa_info['gene_start'],
                'gene_end': aa_info['gene_end'],
                'cds_start': aa_info['cds_start'],
                'cds_end': aa_info['cds_end'],
                'is_silent':aa_info['is_silent'],
                'is_frame_shift': aa_info['is_frame_shift'],
                'ref_aa':aa_info['ref_state'],
                'alt_aa': aa_info['alt_state'],
                'variant_start': non_gap_start,
                'variant_end': non_gap_end,
                'kmer_start': start + 1,
                'kmer_end': end + 1,
                'variant_pos': variant_pos,
                'variant_neg': variant_neg,
                'positive': ''.join(pos_kmer),
                'negative': ''.join(neg_kmer),
                'positive_seqs': [],
                'partial_positive_seqs': []}
        scheme[kmer_name]['positive_seqs'].append(seq_id)

    return scheme


def indel_kmer(input_alignment,indels,consensus_seq,reference_info,ref_name,min_len,max_len,max_ambig,min_members,min_complexity=0.6,n_threads=1):
    scheme = {}
    ref_non_gap_lookup = generate_non_gap_position_lookup(input_alignment[ref_name])
    ref_len = len(input_alignment[ref_name])
    res = []
    pool = Pool(processes=n_threads)

    for seq_id in indels:
        res.append(pool.apply_async(add_indel_kmer_to_scheme, (indels[seq_id], ref_non_gap_lookup, ref_len, input_alignment,
                                           consensus_seq, reference_info, ref_name,
                               min_len, max_len, max_ambig, min_members, min_complexity, 1)))

    pool.close()
    pool.join()
    for x in res:
        if x._success:
            scheme.update(x.get())
    return scheme

def find_indel_kmers(input_alignment,indels,consensus_seq,reference_info,ref_name,min_len,max_len,max_ambig,min_members,min_complexity=0.6,n_threads=1):
    scheme = {}
    ref_non_gap_lookup = generate_non_gap_position_lookup(input_alignment[ref_name])
    for seq_id in indels:
        for indel in indels[seq_id]:
            (vstart,vend) = indel.split(':')
            vstart = int(vstart)
            vend = int(vend)
            vlen = (vend - vstart)

           #determine variant type
            ref_state = input_alignment[ref_name][vstart:vend +1].replace('-','')
            alt_state = input_alignment[seq_id][vstart:vend +1].replace('-','')

            if ref_state == alt_state:
                continue

            if len(ref_state) > len(alt_state):
                type = 'del'
            else:
                type = 'ins'

            if type == 'ins':
                (start, end) = optimize_kmer(pos=vend,
                                             aln_seqs=input_alignment,
                                             reference_sequence=consensus_seq,
                                             min_length=min_len,
                                             max_length=max_len,
                                             min_members=min_members,
                                             max_ambig=max_ambig,
                                             min_complexity=min_complexity,
                                             n_threads=n_threads)


            else:
                #(start, end) = optimize_kmer(vend, input_alignment[seq_id], min_len, max_len, max_ambig, min_complexity)
                (start, end) = optimize_kmer(pos=vend,
                                             aln_seqs=input_alignment,
                                             reference_sequence=input_alignment[seq_id],
                                             min_length=min_len,
                                             max_length=max_len,
                                             min_members=min_members,
                                             max_ambig=max_ambig,
                                             min_complexity=min_complexity)

            if (start == -1 or end == -1) :
                continue

            if vstart < start:
                if type == 'ins':
                    #params = optimize_kmer(vstart, consensus_seq, min_len, max_len, max_ambig, min_complexity=0.5)

                    params = optimize_kmer(pos=vstart,
                                                 aln_seqs=input_alignment,
                                                 reference_sequence=consensus_seq,
                                                 min_length=min_len,
                                                 max_length=max_len,
                                                 min_members=min_members,
                                                 max_ambig=max_ambig,
                                                 min_complexity=min_complexity)

                else:
                    #params = optimize_kmer(vstart, input_alignment[seq_id], min_len, max_len, max_ambig,
                    #                             min_complexity=0.5)

                    params = optimize_kmer(pos=vstart,
                                                 aln_seqs=input_alignment,
                                                 reference_sequence=input_alignment[seq_id],
                                                 min_length=min_len,
                                                 max_length=max_len,
                                                 min_members=min_members,
                                                 max_ambig=max_ambig,
                                                 min_complexity=min_complexity)
                start = params[0]

            rel_start = (vstart - start)
            rel_end = rel_start + vlen

            if rel_end > end:
                rel_end = params[1]

            if (start == -1 or end == -1) :
                continue

            neg_kmer = list(consensus_seq[start:end+1])
            pos_kmer = list(consensus_seq[start:end+1])

            variant_pos = []
            variant_neg = []

            if type == 'ins':
                for i in range(rel_start,rel_end):

                    neg_kmer[i] = '-'
                    variant_neg.append('-')
                    variant_pos.append(pos_kmer[i])

            if type == 'del':
                for i in range(rel_start,rel_end):

                    pos_kmer[i] = '-'
                    variant_neg.append(neg_kmer[i])
                    variant_pos.append('-')

            variant_pos = ''.join(variant_pos)
            variant_neg = ''.join(variant_neg)

            #Handle long indels which exceed the max length of the kmer
            pos_kmer = ''.join(pos_kmer)
            neg_kmer = ''.join(neg_kmer)
            neg_kmer = neg_kmer.replace('-','')
            pos_kmer = pos_kmer.replace('-', '')
            if len(neg_kmer) > max_len:
                diff = False
                for i in range(0,len(neg_kmer)):
                    if i >= len(pos_kmer):
                        break
                    if neg_kmer[i] != pos_kmer[i]:
                        diff = True
                        end = start + i
                    if diff and i >= max_len:
                        break
                if (end - start) < min_len:
                    end = start + min_len
                    i = min_len

            neg_kmer = neg_kmer[0:i+1]
            pos_kmer = pos_kmer[0:i+1]

            if len(neg_kmer) > min_len and len(pos_kmer) > min_len:
                #Trim any degenerate sites at the end of the sequence as long as that isn't a difference
                if neg_kmer[0] == pos_kmer[0]:
                    if neg_kmer[0] not in ['A','T','C','G']:
                        start-=1
                        neg_kmer = neg_kmer[1:len(neg_kmer)]
                        pos_kmer = pos_kmer[1:len(neg_kmer)]
                if len(neg_kmer) == len(pos_kmer):
                    if neg_kmer[len(neg_kmer)-1] == pos_kmer[len(pos_kmer)-1]:
                        if neg_kmer[len(neg_kmer)-1] not in ['A','T','C','G']:
                            end -=1
                            neg_kmer = neg_kmer[0:len(neg_kmer)-1]
                            pos_kmer = pos_kmer[0:len(pos_kmer)-1]

            neg_kmer = list(neg_kmer)
            pos_kmer = list(pos_kmer)
            if variant_pos == variant_neg:
                continue

            non_gap_start = ref_non_gap_lookup[vstart]
            pos = vstart
            while non_gap_start == -1:
                pos -= 1
                non_gap_start = ref_non_gap_lookup[pos]

            non_gap_end  = non_gap_start + vlen


            kmer_name = "{}{}_{}".format(type,non_gap_start, non_gap_end)
            aa_info = get_aa_delta(non_gap_start, non_gap_end, variant_pos, reference_info,ref_name,trans_table=1)

            if len(aa_info['gene']) != 0:
                kmer_name_aa = "{}{}{}".format(aa_info['ref_state'], aa_info['aa_start'], aa_info['alt_state'])
            else:
                kmer_name_aa = 'N/A'
            if kmer_name not in scheme:
                scheme[kmer_name] = {
                    'dna_name': kmer_name,
                    'aa_name': kmer_name_aa,
                    'type': type,
                    'gene':aa_info['gene'],
                    'gene_start':aa_info['gene_start'],
                    'gene_end': aa_info['gene_end'],
                    'cds_start': aa_info['cds_start'],
                    'cds_end': aa_info['cds_end'],
                    'is_silent':aa_info['is_silent'],
                    'is_frame_shift': aa_info['is_frame_shift'],
                    'ref_aa':aa_info['ref_state'],
                    'alt_aa': aa_info['alt_state'],
                    'variant_start': non_gap_start,
                    'variant_end': non_gap_end,
                    'kmer_start': start + 1,
                    'kmer_end': end + 1,
                    'variant_pos': variant_pos,
                    'variant_neg': variant_neg,
                    'positive': ''.join(pos_kmer),
                    'negative': ''.join(neg_kmer),
                    'positive_seqs': [],
                    'partial_positive_seqs': []}
            scheme[kmer_name]['positive_seqs'].append(seq_id)

    return scheme



def find_indel_kmers_bck(input_alignment,indels,consensus_seq,reference_info,ref_name,min_len,max_len,max_ambig,min_members,min_complexity=0.6,n_threads=1):
    scheme = {}
    ref_non_gap_lookup = generate_non_gap_position_lookup(input_alignment[ref_name])
    for seq_id in indels:
        for indel in indels[seq_id]:
            (vstart,vend) = indel.split(':')
            vstart = int(vstart)
            vend = int(vend)
            vlen = (vend - vstart)

           #determine variant type
            ref_state = input_alignment[ref_name][vstart:vend +1].replace('-','')
            alt_state = input_alignment[seq_id][vstart:vend +1].replace('-','')

            if ref_state == alt_state:
                continue

            if len(ref_state) > len(alt_state):
                type = 'del'
            else:
                type = 'ins'

            if type == 'ins':
                (start, end) = optimize_kmer(pos=vend,
                                             aln_seqs=input_alignment,
                                             reference_sequence=consensus_seq,
                                             min_length=min_len,
                                             max_length=max_len,
                                             min_members=min_members,
                                             max_ambig=max_ambig,
                                             min_complexity=min_complexity,
                                             n_threads=n_threads)


            else:
                #(start, end) = optimize_kmer(vend, input_alignment[seq_id], min_len, max_len, max_ambig, min_complexity)
                (start, end) = optimize_kmer(pos=vend,
                                             aln_seqs=input_alignment,
                                             reference_sequence=input_alignment[seq_id],
                                             min_length=min_len,
                                             max_length=max_len,
                                             min_members=min_members,
                                             max_ambig=max_ambig,
                                             min_complexity=min_complexity)

            if (start == -1 or end == -1) :
                continue

            if vstart < start:
                if type == 'ins':
                    #params = optimize_kmer(vstart, consensus_seq, min_len, max_len, max_ambig, min_complexity=0.5)

                    params = optimize_kmer(pos=vstart,
                                                 aln_seqs=input_alignment,
                                                 reference_sequence=consensus_seq,
                                                 min_length=min_len,
                                                 max_length=max_len,
                                                 min_members=min_members,
                                                 max_ambig=max_ambig,
                                                 min_complexity=min_complexity)

                else:
                    #params = optimize_kmer(vstart, input_alignment[seq_id], min_len, max_len, max_ambig,
                    #                             min_complexity=0.5)

                    params = optimize_kmer(pos=vstart,
                                                 aln_seqs=input_alignment,
                                                 reference_sequence=input_alignment[seq_id],
                                                 min_length=min_len,
                                                 max_length=max_len,
                                                 min_members=min_members,
                                                 max_ambig=max_ambig,
                                                 min_complexity=min_complexity)
                start = params[0]

            rel_start = (vstart - start)
            rel_end = rel_start + vlen

            if rel_end > end:
                rel_end = params[1]

            if (start == -1 or end == -1) :
                continue

            neg_kmer = list(consensus_seq[start:end+1])
            pos_kmer = list(consensus_seq[start:end+1])

            variant_pos = []
            variant_neg = []

            if type == 'ins':
                for i in range(rel_start,rel_end):

                    neg_kmer[i] = '-'
                    variant_neg.append('-')
                    variant_pos.append(pos_kmer[i])

            if type == 'del':
                for i in range(rel_start,rel_end):

                    pos_kmer[i] = '-'
                    variant_neg.append(neg_kmer[i])
                    variant_pos.append('-')

            variant_pos = ''.join(variant_pos)
            variant_neg = ''.join(variant_neg)

            #Handle long indels which exceed the max length of the kmer
            pos_kmer = ''.join(pos_kmer)
            neg_kmer = ''.join(neg_kmer)
            neg_kmer = neg_kmer.replace('-','')
            pos_kmer = pos_kmer.replace('-', '')
            if len(neg_kmer) > max_len:
                diff = False
                for i in range(0,len(neg_kmer)):
                    if i >= len(pos_kmer):
                        break
                    if neg_kmer[i] != pos_kmer[i]:
                        diff = True
                        end = start + i
                    if diff and i >= max_len:
                        break
                if (end - start) < min_len:
                    end = start + min_len
                    i = min_len

            neg_kmer = neg_kmer[0:i+1]
            pos_kmer = pos_kmer[0:i+1]

            if len(neg_kmer) > min_len and len(pos_kmer) > min_len:
                #Trim any degenerate sites at the end of the sequence as long as that isn't a difference
                if neg_kmer[0] == pos_kmer[0]:
                    if neg_kmer[0] not in ['A','T','C','G']:
                        start-=1
                        neg_kmer = neg_kmer[1:len(neg_kmer)]
                        pos_kmer = pos_kmer[1:len(neg_kmer)]
                if len(neg_kmer) == len(pos_kmer):
                    if neg_kmer[len(neg_kmer)-1] == pos_kmer[len(pos_kmer)-1]:
                        if neg_kmer[len(neg_kmer)-1] not in ['A','T','C','G']:
                            end -=1
                            neg_kmer = neg_kmer[0:len(neg_kmer)-1]
                            pos_kmer = pos_kmer[0:len(pos_kmer)-1]

            neg_kmer = list(neg_kmer)
            pos_kmer = list(pos_kmer)
            if variant_pos == variant_neg:
                continue

            non_gap_start = ref_non_gap_lookup[vstart]
            pos = vstart
            while non_gap_start == -1:
                pos -= 1
                non_gap_start = ref_non_gap_lookup[pos]

            non_gap_end  = non_gap_start + vlen


            kmer_name = "{}{}_{}".format(type,non_gap_start, non_gap_end)
            aa_info = get_aa_delta(non_gap_start, non_gap_end, variant_pos, reference_info,ref_name,trans_table=1)

            if len(aa_info['aa_start']) != 0:
                kmer_name_aa = "{}{}{}".format(aa_info['ref_state'], aa_info['aa_start'], aa_info['alt_state'])
            else:
                kmer_name_aa = 'N/A'
            if kmer_name not in scheme:
                scheme[kmer_name] = {
                    'dna_name': kmer_name,
                    'aa_name': kmer_name_aa,
                    'type': type,
                    'gene':aa_info['gene'],
                    'gene_start':aa_info['gene_start'],
                    'gene_end': aa_info['gene_end'],
                    'cds_start': aa_info['cds_start'],
                    'cds_end': aa_info['cds_end'],
                    'is_silent':aa_info['is_silent'],
                    'is_frame_shift': aa_info['is_frame_shift'],
                    'ref_aa':aa_info['ref_state'],
                    'alt_aa': aa_info['alt_state'],
                    'variant_start': non_gap_start,
                    'variant_end': non_gap_end,
                    'kmer_start': start + 1,
                    'kmer_end': end + 1,
                    'variant_pos': variant_pos,
                    'variant_neg': variant_neg,
                    'positive': ''.join(pos_kmer),
                    'negative': ''.join(neg_kmer),
                    'positive_seqs': [],
                    'partial_positive_seqs': []}
            scheme[kmer_name]['positive_seqs'].append(seq_id)

    return scheme


def count_ambig(seq):
    # count ambiguous characters
    bases = ['A', 'T', 'C', 'G']
    nt_count = 0
    mers = count_kmers(seq, K=1)
    for b in bases:
        if b in mers:
            nt_count += mers[b]
    return (len(seq) - nt_count)


def get_kmer_complexity(kmer,K=2):
    # determine the complexity of the sequence and remove kmers composed heavily of the same 2-mer
    mers = count_kmers(kmer, K)
    num_mers = sum(mers.values())
    mer_perc = []

    for m in mers:
        mer_perc.append(mers[m] / num_mers)

    if len(mer_perc) == 0:
        return {'average':0,'min':0,'max':0}
    else:
        return {'average': sum(mer_perc)/len(mer_perc), 'min': min(mer_perc), 'max': max(mer_perc)}


def get_kmer_groups(scheme):
    groups = {}
    for kmer_id in scheme:
        group = scheme[kmer_id]['group_id']
        if not group in groups:
            groups[group] = []
        groups[group].append(str(scheme[kmer_id]['key']))
    return groups

def get_kmer_group_mapping(scheme):
    groups = {}
    for kmer_id in scheme:
        group = scheme[kmer_id]['group_id']
        groups[str(scheme[kmer_id]['key'])] = group
    return groups

def add_key(scheme):
    i = 0
    for kmer_id in scheme:
        scheme[kmer_id]['key'] = i
        i+=1
    return scheme

def is_seq_in(target,query):
    len_target = len(target)
    len_query = len(query)
    if len_target != len_query :
        return False
    for i in range(0,len_query ):
        if target[i] in bases_dict:
            baseT = bases_dict[target[i]]
        else:
            baseT = [target[i]]
        if query[i] in bases_dict:
            baseQ = bases_dict[query[i]]
        else:
            baseQ = [query[i]]
        if baseQ != baseT:
            intersect = set(baseT) & set(baseQ)
            if len(intersect) == 0:
                return False

    return True



def process_biohansel_kmer(scheme_kmer_groups,scheme_target_to_group_mapping,scheme,kmer_df,min_cov=20):
    data = {}
    columns = kmer_df.columns
    if not 'seq' in columns:
        return data
    for row in kmer_df.itertuples():
        seq = row.seq
        target_id = str(row.refposition)
        is_pos_kmer = row.is_pos_kmer
        sample = row.sample
        if 'freq' in columns:
            freq = row.freq
        else:
            freq = min_cov
        if not sample in data:
            data[sample] = {
            }
        if not target_id in data[sample]:
            data[sample][target_id] = {}

        if not seq in data[sample][target_id]:
            data[sample][target_id][seq] = {
                'freq': freq,
                'is_pos_kmer': is_pos_kmer
            }
        else:
            data[sample][target_id][seq]['freq']+= freq

        if target_id in scheme_target_to_group_mapping:
            group_id = scheme_target_to_group_mapping[target_id]
            if group_id in scheme_kmer_groups:
                members = scheme_kmer_groups[group_id]
                if len(members) == 1:
                    continue
                for member in members:
                    if target_id == member:
                        continue
                    positive = scheme[member]['positive']
                    negative = scheme[member]['negative']
                    is_present = False
                    if is_seq_in(seq, positive):
                        is_present = True
                        is_pos_kmer = True
                    elif is_seq_in(seq, negative):
                        is_present = True
                        is_pos_kmer = False
                    if is_present:
                        if not member in data[sample]:
                            data[sample][member] = {}

                        if not seq in data[sample][member]:
                            data[sample][member][seq] = {
                                'freq': freq,
                                'is_pos_kmer': is_pos_kmer
                            }
                        else:
                            data[sample][member][seq]['freq'] += freq
    return data


def sliding_window_kmer(seq,kmer_size):
    kmers = []
    l = len(seq) - kmer_size
    for i in range(0,l):
        kmers.append(seq[i:i+kmer_size])
    return kmers

def generate_ref_kmers(aln,min_ksize,max_k_size):
    seqs = {}
    for seq_id in aln:
        if not seq_id in seqs:
            seqs[seq_id] = {}
        seq = aln[seq_id].replace('-','')
        for i in range(min_ksize,max_k_size):
            if not i in seqs[seq_id]:
                seqs[seq_id][i] = []
            seqs[seq_id][i].extend(sliding_window_kmer(seq, i))
            seqs[seq_id][i].extend(sliding_window_kmer(revcomp(seq), i))
    return seqs

def get_count_kmer_occurances(search,sample_kmers):
    counts = {}
    for sample in sample_kmers:
        counts[sample] = 0
        for k in sample_kmers[sample]:
            if k == search:
                counts[sample]+=1
    return counts

def is_kmer_valid(df,min_members):
    is_valid = True
    counts = df['contig_id'].describe()
    if counts['freq'] > 1:
        is_valid = False
    return is_valid

def read_genotype_profiles(file):
    data = read_tsv(file)
    profiles = {}
    for row in data.itertuples():
        sample_id = row.sample_id
        genotype = row.genotype
        profile = set(row.kmer_profile.replace(' ','').split(','))
        md5 = row.kmer_st_md5
        phrase = row.kmer_phrase
        profiles[sample_id] = {
            'genotype':genotype,
            'profile':profile,
            'md5':md5,
            'phrase':phrase
        }
    return profiles

def dist_compatible_profiles(profiles):
    data = []
    for sample_id in profiles:
        data.append(profiles[sample_id]['profile'])
    return data

def generate_target_presence_table(sample_kmer_data,min_ratio,max_ratio):
    profiles = {}
    for sample_id in sample_kmer_data:
        profiles[sample_id] = {}
        data = sample_kmer_data[sample_id]['ratios']
        for target in data:

            value = data[target]
            if value == -1:
                value = float("NaN")
            elif value >= max_ratio:
                value = 1
            elif value >= min_ratio:
                value = 0.5
            else:
                value = 0
            profiles[sample_id][target] = value
    return pd.DataFrame.from_dict(profiles,orient='columns')

def nan_compatible_kmer_pairwise_distmatrix(profile_df):
    '''
    Computes pairwise jaccard distances between sample profiles
    '''
    i = profile_df.values.T
    j = np.nansum((i - i[:, None]) ** 2, axis=2) ** .5
    #df = (lambda v, c: pd.DataFrame(v, c, c))(j, profile_df.columns)
    return j



def summarize_kmer_targets(scheme_kmer_target_keys,sample_kmer_data,min_cov):
    kmer_data = {}
    for kmer_id in scheme_kmer_target_keys:
        kmer_data[kmer_id] = {'count_negative': [], 'count_positive': [], 'missing_samples': [], 'count_found': 0}
    for sample_id in sample_kmer_data:
        counts = sample_kmer_data[sample_id]['counts']
        for kmer_id in counts:
            pos_count = counts[kmer_id]['positive']
            neg_count = counts[kmer_id]['negative']
            total_count = pos_count + neg_count
            if total_count >= min_cov:
                kmer_data[kmer_id]['count_found'] += 1
                kmer_data[kmer_id]['count_positive'].append(pos_count)
                kmer_data[kmer_id]['count_negative'].append(neg_count)
            else:
                kmer_data[kmer_id]['missing_samples'].append(sample_id)
    return kmer_data

def get_valid_targets(kmer_data,num_samples,min_pos_freq,max_frac_missing):
    valid = []
    for target in kmer_data:
        pos_count = len(kmer_data[target]['count_positive'])
        found_count = kmer_data[target]['count_found']
        if found_count > 0:
            missing_frac = 1 - (found_count/num_samples)
        status = 'Pass'
        if pos_count < min_pos_freq:
            status = 'Fail'
        elif missing_frac > max_frac_missing:
            status = 'Fail'
        if status == 'Pass':
            valid.append(target)
    return valid

def get_kmer_freq_by_genotype(sample_kmer_data,reported_genotypes,min_cov,min_cov_frac):
    genotype_data = {}
    for sample_id in sample_kmer_data:
        counts = sample_kmer_data[sample_id]['counts']
        ratios = sample_kmer_data[sample_id]['ratios']
        genotype = reported_genotypes[sample_id]
        if not genotype in genotype_data:
            genotype_data[genotype] = {}
        for kmer_id in counts:
            if kmer_id not in genotype_data[genotype]:
                genotype_data[genotype][kmer_id] = {
                    'positive': [],
                    'negative': [],
                    'missing': []
                }
            pos_count = counts[kmer_id]['positive']
            neg_count = counts[kmer_id]['negative']
            total_count = pos_count + neg_count
            if total_count < min_cov:
                genotype_data[genotype][kmer_id]['missing'].append(sample_id)
            else:
                if pos_count >= min_cov and ratios[kmer_id] >= min_cov_frac:
                    genotype_data[genotype][kmer_id]['positive'].append(sample_id)
                if neg_count >= min_cov and 1-ratios[kmer_id] >= min_cov_frac:
                    genotype_data[genotype][kmer_id]['negative'].append(sample_id)


    return genotype_data

def process_rules(sample_kmer_data,valid_targets,rules,genotypes,reported_genotypes,min_cov,min_cov_frac):
    scores = {}
    genotypes = set(genotypes)

    for sample_id in sample_kmer_data:
        counts = sample_kmer_data[sample_id]['counts']
        ratios = sample_kmer_data[sample_id]['ratios']
        genotype = reported_genotypes[sample_id]
        scores[sample_id] = {
            'conflicting_targets':[],
            'reported_genotype':'',
            'candidates':[],
            'exclude':[],
            'include':[]
        }
        for kmer_id in counts:
            if not kmer_id in rules or len(rules[kmer_id]['positive']) == 0 or kmer_id not in valid_targets:
                continue
            pos_count = counts[kmer_id]['positive']
            neg_count = counts[kmer_id]['negative']
            total_count = pos_count + neg_count
            scores[sample_id]['reported_genotype'] = reported_genotypes[sample_id]
            if total_count < min_cov:
                continue
            is_mixed = pos_count >= min_cov and neg_count > min_cov and ratios[kmer_id] >= min_cov_frac
            if is_mixed:
                continue
            if pos_count >= min_cov and ratios[kmer_id] >= min_cov_frac:
                exclude = list(genotypes - (set(rules[kmer_id]['positive']).update(set(rules[kmer_id]['partials']))))
                scores[sample_id]['exclude'].append(exclude)

            else:
                exclude = list(set(rules[kmer_id]['positive']))
                scores[sample_id]['exclude'].append(exclude)

            if genotype in exclude:
                scores[sample_id]['conflicting_targets'].append(kmer_id)
        scores[sample_id]['candidates'] = list(genotypes - set(scores[sample_id]['exclude']))
    return scores

def evaluate_rules(scores,rules,threshold):
    results = {
        'scheme_score':0,
        'genotype_scores':{},
        'rules':rules
    }
    conflicting_kmers = {}
    genotype_counts = {}
    genotype_scores = {}

    for sample_id in scores:
        scores[sample_id]['conflicting_targets'] = list(set(scores[sample_id]['conflicting_targets']))
        genotype = scores[sample_id]['reported_genotype']
        if not genotype in genotype_counts:
            genotype_counts[genotype] = 0
            genotype_scores[genotype] = 0
        genotype_counts[genotype] += 1

        if genotype in scores[sample_id]['candidates']:
            genotype_scores[genotype] += 1 / len(scores[sample_id]['candidates'])

        for kmer_id in scores[sample_id]['conflicting_targets']:
            if not kmer_id in conflicting_kmers:
                conflicting_kmers[kmer_id] = {}
            if not genotype in conflicting_kmers[kmer_id]:
                conflicting_kmers[kmer_id][genotype] = 0
            conflicting_kmers[kmer_id][genotype] +=1
    scheme_score = 0
    for genotype in genotype_scores:
        scheme_score+= genotype_scores[genotype]
    num_genotypes = len(genotype_scores)
    results['scheme_score'] = scheme_score
    results['genotype_scores'] = genotype_scores
    num_samples = len(scores)
    for kmer_id in conflicting_kmers:
        kRules = rules[kmer_id]
        sum_conflicts = 0
        for genotype in conflicting_kmers[kmer_id]:
            sum_conflicts += conflicting_kmers[kmer_id][genotype]
            percent_conflicting = conflicting_kmers[kmer_id][genotype] / genotype_counts[genotype]
            if percent_conflicting < threshold:
                continue
            if genotype in kRules['positive']:
                kRules['positive'] = list(set(kRules['positive']) - set([genotype]))
            kRules['partials'].append(genotype)

        if sum_conflicts / num_samples > threshold or \
            len(kRules['positive']) == len(num_genotypes):

            kRules['positive'] = []
            kRules['partials'] = []

        rules[kmer_id] = kRules
    results['rules'] = rules
    return results


