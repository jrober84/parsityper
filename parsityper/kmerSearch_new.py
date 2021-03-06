from Bio.Seq import Seq
from Bio import SeqIO
import itertools
from collections import defaultdict
import sys,re
import pandas as pd
from ahocorasick import Automaton
from itertools import product
from typing import List, Any, Optional, Tuple, Union
NT_SUB = str.maketrans('acgtrymkswhbvdnxACGTRYMKSWHBVDNX',
                       'tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX')
REGEX_GZIPPED = re.compile(r'^.+\.gz$')
import time
import re
from multiprocessing import Pool

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


def init_automaton(scheme_fasta):
    """Initialize Aho-Corasick Automaton with kmers from SNV scheme fasta

    Args:
        scheme_fasta: SNV scheme fasta file path

    Returns:
         Aho-Corasick Automaton with kmers loaded
    """
    A = Automaton()
    for header, sequence in parse_fasta(scheme_fasta):
        kmer_list = expand_degenerate_bases(sequence)
        for idx,seq in enumerate(kmer_list):
            A.add_word(seq, (header, seq, False))
            A.add_word(revcomp(seq), (header, seq, True))
    A.make_automaton()
    return A

def init_automaton_dict(seqs):
    """Initialize Aho-Corasick Automaton with kmers from SNV scheme fasta

    Args:
        scheme_fasta: SNV scheme fasta file path

    Returns:
         Aho-Corasick Automaton with kmers loaded
    """
    stime = time.time()
    A = Automaton()
    for seq_id in seqs:
        sequence = seqs[seq_id]
        kmer_list = expand_degenerate_bases(sequence)
        for idx,seq in enumerate(kmer_list):
            A.add_word(seq, (seq_id, seq, False))
            A.add_word(revcomp(seq), (seq_id, seq, True))
    A.make_automaton()
    print("time taken for build: {}".format(time.time() - stime))
    return A

def parse_fasta(filepath):
    """Parse a FASTA/FASTA.GZ file returning a generator yielding tuples of fasta headers to sequences.

    Args:
        filepath (str): Fasta file path

    Returns:
        generator: yields tuples of (<fasta header>, <fasta sequence>)
    """
    if REGEX_GZIPPED.match(filepath):
        # assumes Linux os with zcat installed
        import os
        with os.popen('zcat < {}'.format(filepath)) as f:
            yield from SimpleFastaParser(f)
    else:
        with open(filepath, 'r') as f:
            yield from SimpleFastaParser(f)

def find_in_fasta(automaton: Automaton, fasta: str) -> pd.DataFrame:
    """Find scheme kmers in input fasta file

    Args:
        automaton: Aho-Corasick Automaton with scheme SNV target kmers loaded
        fasta: Input fasta path

    Returns:
        Dataframe with any matches found in input fasta file
    """
    res = []
    for contig_header, sequence in parse_fasta(fasta):
        for idx, (kmername, kmer_seq, is_revcomp) in automaton.iter(sequence):
            res.append((kmername, kmer_seq, is_revcomp, contig_header, idx))
    columns = ['kmername', 'seq', 'is_revcomp', 'contig_id', 'match_index']
    return pd.DataFrame(res, columns=columns)

def find_in_fasta_dict(automaton: Automaton, seqs: dict) -> pd.DataFrame:
    """Find scheme kmers in input fasta file

    Args:
        automaton: Aho-Corasick Automaton with scheme SNV target kmers loaded
        fasta: Input fasta path

    Returns:
        Dataframe with any matches found in input fasta file
    """
    res = []
    for seq_id in seqs:
        seq = seqs[seq_id].replace('-','')
        for idx, (kmername, kmer_seq, is_revcomp) in automaton.iter(seq):
            res.append((kmername, kmer_seq, is_revcomp, seq_id, idx))
    columns = ['kmername', 'seq', 'is_revcomp', 'contig_id', 'match_index']
    return pd.DataFrame(res, columns=columns)


def parallel_query_contigs(input_genomes,
                            automaton: Automaton,
                           n_threads: int = 1):
    print("===>{}".format(n_threads))
    stime = time.time()
    if n_threads == 1:
        return find_in_fasta_dict(automaton,input_genomes)
    else:

        pool = Pool(processes=n_threads)
        res = []
        it = iter(input_genomes)
        length = len(input_genomes)
        chunk_size = int(length / n_threads)

        for i in range(0, length, chunk_size):
            chunk = {}
            sub = time.time()
            for seq_id in itertools.islice(it,chunk_size):
                seq = input_genomes[seq_id]
                chunk[seq_id] = seq
            res.append(pool.apply_async(find_in_fasta_dict,  ( automaton, chunk  )))
            print("time taken to submit: {}".format(time.time() - sub))
        #cleanup
        pool.close()
        pool.join()
        pd.concat([x.get() for x in res])
        print("time taken for search: {} with {} threads".format(time.time()-stime,n_threads))
        stime = time.time()
        find_in_fasta_dict(automaton, input_genomes)
        print("time taken for search: {} with {} threads".format(time.time() - stime, 1))

        return pd.concat([x.get() for x in res])


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


def find_in_fastqs(automaton: Automaton, *fastqs):
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
    for kmer_seq, freq in kmer_seq_counts.items():
        kmername, sequence, _ = automaton.get(kmer_seq)
        res.append((kmername, kmer_seq, freq))
    return pd.DataFrame(res, columns=['kmername', 'seq', 'freq'])
