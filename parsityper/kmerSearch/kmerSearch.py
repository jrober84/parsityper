import itertools
import os.path
from collections import defaultdict
import sys
import pandas as pd
from ahocorasick import Automaton
from itertools import product
from typing import List, Any, Optional, Tuple, Union
import time
import re
from multiprocessing import Pool, current_process
from os import fork, getpid

bases_dict = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T'], }

REGEX_GZIPPED = re.compile(r'^.+\.gz$')

NT_SUB = str.maketrans('acgtrymkswhbvdnxACGTRYMKSWHBVDNX',
                       'tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX')

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


def init_automaton_dict(seqs):
    """Initialize Aho-Corasick Automaton with kmers from SNV scheme fasta

    Args:
        scheme_fasta: SNV scheme fasta file path

    Returns:
         Aho-Corasick Automaton with kmers loaded
    """
    A = Automaton()
    for seq_id in seqs:
        sequence = seqs[seq_id]
        kmer_list = expand_degenerate_bases(sequence.replace('-',''))
        for idx,seq in enumerate(kmer_list):
            A.add_word(seq, (seq_id, seq, False))
            A.add_word(revcomp(seq), (seq_id, seq, True))
    A.make_automaton()
    return A

def find_in_fasta_dict(automaton: Automaton, seqs: dict) -> pd.DataFrame:
    """Find scheme kmers in input fasta file

    Args:
        automaton: Aho-Corasick Automaton with scheme SNV target kmers loaded
        fasta: Input fasta path

    Returns:
        Dataframe with any matches found in input fasta file
    """
    res = []
    iter_keys = seqs.keys()
    for seq_id in iter_keys:
        seq = seqs[seq_id].replace('-','')
        for idx, (kmername, kmer_seq, is_revcomp) in automaton.iter(seq):
            res.append((kmername, kmer_seq, is_revcomp, seq_id, idx))
    columns = ['kmername', 'seq', 'is_revcomp', 'contig_id', 'match_index']
    return pd.DataFrame(res, columns=columns)

def parse_fastq(filepath):
    """Parse a FASTQ/FASTQ.GZ file returning a generator yielding tuples of FASTQ entry headers and sequences.

    Args:
        filepath (str): FASTQ/FASTQ.GZ file path

    Returns:
        generator: yields tuples of (<fastq header>, <fastq sequence>)
    """
    if REGEX_GZIPPED.match(filepath):
        # using os.popen with zcat since it is much faster than gzip.open or gzip.open(io.BufferedReader)
        # http://aripollak.com/pythongzipbenchmarks/
        # assumes Linux os with zcat installed
        import os
        with os.popen('zcat < {}'.format(filepath)) as f:
            yield from _parse_fastq(f)
    else:
        with open(filepath, 'r') as f:
            yield from _parse_fastq(f)

def _parse_fastq(f):
    """Simple FASTQ parser which yields the header and sequence ignoring the quality scores

    Args:
        f: file-like object

    Yields:
        Tuple of FASTQ entry header and sequence
    """
    header = ''
    seq = ''
    skip = False
    for line in f:
        if skip:
            skip = False
            continue
        line = line.strip()
        if line == '':
            continue
        if line[0] == '@':
            header = line.replace('@', '')
        elif line[0] == '+':
            yield header, seq
            skip = True
        else:
            seq = line.upper()

def find_in_fastqs(automaton: Automaton, fastqs):
    """Find scheme kmers in input fastq files

    Args:
        automaton: Aho-Corasick Automaton with scheme SNV target kmers loaded
        fastqs: Input fastq file paths

    Returns:
        Dataframe with any matches found in input fastq files
    """
    kmer_seq_counts = defaultdict(int)
    for fastq in fastqs:
        for _, sequence in parse_fastq(fastq):
            for idx, (_, kmer_seq, _) in automaton.iter(sequence):
                kmer_seq_counts[kmer_seq] += 1
    res = []
    kmer_freq = {}
    for kmer_seq, freq in kmer_seq_counts.items():
        kmername, sequence, _ = automaton.get(kmer_seq)
        if not kmername in kmer_freq:
            kmer_freq[kmername] = 0
        kmer_freq[kmername]+= freq
        res.append((kmername, kmer_seq, freq))

    return pd.DataFrame(res, columns=['kmername','kmer_seq','freq'])

def perform_kmerSearch_fastq(uid_to_kseq_map,kseq_to_uids,aho,sequence_files):
    kmer_df = find_in_fastqs(aho, sequence_files)
    kmer_results = {}
    for uid in uid_to_kseq_map:
        kmer_results[uid] = 0
    for row in kmer_df.itertuples():
        uid = row.kmername
        freq = row.freq
        kseq = uid_to_kseq_map[uid]
        uids = kseq_to_uids[kseq]
        for i in uids:
            kmer_results[i] = freq
    return kmer_results

def perform_kmerSearch_fasta(uid_to_kseq_map,kseq_to_uids,aho,sequence_dict,min_freq):
    kmer_df = find_in_fasta_dict(aho, sequence_dict)
    kmer_results = {}
    for uid in uid_to_kseq_map:
        kmer_results[uid] = 0
    for row in kmer_df.itertuples():
        uid = row.kmername
        kseq = uid_to_kseq_map[uid]
        uids = kseq_to_uids[kseq]
        for i in uids:
            kmer_results[i] += min_freq
    return kmer_results

def process_kmer_results(scheme_info,kmer_results,min_freq,min_cov_frac):
    mutation_to_uid = scheme_info['mutation_to_uid']
    uid_to_state = scheme_info['uid_to_state']
    sample_kmer_results = {}
    for sample_id in kmer_results:
        sample_kmer_results[sample_id] = {
            'raw_kmer_freq':kmer_results[sample_id],
            'total_mutation_key_freq':{},
            'positive_mutation_key_freq': {},
            'mutation_frac':{},
            'missing_sites':[],
            'num_missing_sites':0,
            'mixed_sites':[],
            'num_detected_kmers':0,
            'num_mixed_sites': 0,
            'positive_kmers':[],
            'detected_scheme_kmers':[],
            'valid_uids':[],
            'num_positive_kmers': 0,
            'average_kmer_freq':0
        }
        for uid in kmer_results[sample_id]:
            if kmer_results[sample_id][uid] < min_freq:
                kmer_results[sample_id][uid] = 0
    for sample_id in kmer_results:
        kmer_freqs = []
        for mutation_key in mutation_to_uid:
            sample_kmer_results[sample_id]['total_mutation_key_freq'][mutation_key] = 0
            sample_kmer_results[sample_id]['positive_mutation_key_freq'][mutation_key] = 0
            sample_kmer_results[sample_id]['mutation_frac'][mutation_key] = 0
            kmer_uids = mutation_to_uid[mutation_key]

            #multiple ref kmers can map to the same event so take the max value one
            ref_kmer_freq = []
            alt_kmer_freq = []
            ref_kmers_uids = []
            alt_kmers_uids = []
            for uid in kmer_uids:
                freq = kmer_results[sample_id][uid]
                state = uid_to_state[uid]
                if freq < min_freq:
                    kmer_results[sample_id][uid] = 0
                    continue
                if state == 'ref':
                    ref_kmer_freq.append(freq)
                    ref_kmers_uids.append(uid)
                else:
                    alt_kmer_freq.append(freq)
                    alt_kmers_uids.append(uid)
                kmer_freqs.append(freq)

            ref_mutation_freq = 0
            alt_mutation_freq = 0

            if len(set(ref_kmer_freq)) == 1:
                ref_mutation_freq = ref_kmer_freq[0]
            elif len(ref_kmer_freq) > 1:
                ref_mutation_freq = sum(ref_kmer_freq) / len(ref_kmer_freq)

            if len(set(alt_kmer_freq)) == 1:
                alt_mutation_freq = alt_kmer_freq[0]
            elif len(alt_kmer_freq) > 1:
                alt_mutation_freq = sum(alt_kmer_freq) / len(alt_kmer_freq)

            sample_kmer_results[sample_id]['total_mutation_key_freq'][mutation_key] = ref_mutation_freq + alt_mutation_freq
            sample_kmer_results[sample_id]['positive_mutation_key_freq'][mutation_key] = alt_mutation_freq

            if sample_kmer_results[sample_id]['total_mutation_key_freq'][mutation_key] > 0:
                sample_kmer_results[sample_id]['mutation_frac'][mutation_key] = alt_mutation_freq / (ref_mutation_freq + alt_mutation_freq)


            if  sample_kmer_results[sample_id]['total_mutation_key_freq'][mutation_key] > 0:
                sample_kmer_results[sample_id]['average_kmer_freq'] = sum(kmer_freqs) / len(kmer_freqs)

            freq = sample_kmer_results[sample_id]['total_mutation_key_freq'][mutation_key]
            if freq < min_freq:
                sample_kmer_results[sample_id]['missing_sites'].append(mutation_key)

            frac = sample_kmer_results[sample_id]['mutation_frac'][mutation_key]

            if freq >= min_freq and frac >= min_cov_frac and frac <= 1 - min_cov_frac:
                sample_kmer_results[sample_id]['mixed_sites'].append(mutation_key)
                sample_kmer_results[sample_id]['detected_scheme_kmers'].extend(kmer_uids)
            elif freq >= min_freq:
                if frac >= min_cov_frac:
                    sample_kmer_results[sample_id]['detected_scheme_kmers'].extend(alt_kmers_uids)
                    for uid in ref_kmers_uids:
                        kmer_results[sample_id][uid] = 0

                else:
                    sample_kmer_results[sample_id]['detected_scheme_kmers'].extend(ref_kmers_uids)
                    for uid in alt_kmers_uids:
                        kmer_results[sample_id][uid] = 0

            pos_freq = sample_kmer_results[sample_id]['positive_mutation_key_freq'][mutation_key]
            if pos_freq > min_freq and frac >= min_cov_frac:
                kmer_uids = mutation_to_uid[mutation_key]
                for uid in kmer_uids:
                    f = kmer_results[sample_id][uid]
                    state = uid_to_state[uid]
                    if state == 'alt':
                        if f >= min_freq:
                            sample_kmer_results[sample_id]['positive_kmers'].append(uid)
        for mutation_key in mutation_to_uid:
            if mutation_key in sample_kmer_results[sample_id]['missing_sites']:
                continue
            sample_kmer_results[sample_id]['valid_uids'].extend(mutation_to_uid[mutation_key])
        sample_kmer_results[sample_id]['valid_uids'] = list(set(sample_kmer_results[sample_id]['valid_uids']))
        sample_kmer_results[sample_id]['num_positive_kmers'] = len(sample_kmer_results[sample_id]['positive_kmers'])
        sample_kmer_results[sample_id]['num_missing_sites'] = len(sample_kmer_results[sample_id]['missing_sites'])
        sample_kmer_results[sample_id]['num_mixed_sites'] = len(sample_kmer_results[sample_id]['mixed_sites'])
        sample_kmer_results[sample_id]['num_detected_kmers'] = len(kmer_freqs)
        sample_kmer_results[sample_id]['detected_scheme_kmers'] = list(set(sample_kmer_results[sample_id]['detected_scheme_kmers']))
    return sample_kmer_results
