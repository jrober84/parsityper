import itertools
import os.path
from collections import defaultdict
import sys
import pandas as pd
from ahocorasick import Automaton
from itertools import product
from typing import List, Any, Optional, Tuple, Union
NT_SUB = str.maketrans('acgtrymkswhbvdnxACGTRYMKSWHBVDNX',
                       'tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX')
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
    return A


def SimpleFastaParser(handle):
    """Iterate over Fasta records as string tuples.
    Arguments:
     - handle - input stream opened in text mode
    For each record a tuple of two strings is returned, the FASTA title
    line (without the leading '>' character), and the sequence (with any
    whitespace removed). The title line is not divided up into an
    identifier (the first word) and comment or description.
    >>> with open("Fasta/dups.fasta") as handle:
    ...     for values in SimpleFastaParser(handle):
    ...         print(values)
    ...
    ('alpha', 'ACGTA')
    ('beta', 'CGTC')
    ('gamma', 'CCGCC')
    ('alpha (again - this is a duplicate entry to test the indexing code)', 'ACGTA')
    ('delta', 'CGCGC')
    """
    # Skip any text before the first record (e.g. blank lines, comments)
    for line in handle:
        if line[0] == ">":
            title = line[1:].rstrip()
            break
    else:
        # no break encountered - probably an empty file
        return

    # Main logic
    # Note, remove trailing whitespace, and any internal spaces
    # (and any embedded \r which are possible in mangled files
    # when not opened in universal read lines mode)
    lines = []
    for line in handle:
        if line[0] == ">":
            yield title, "".join(lines).replace(" ", "").replace("\r", "").upper()
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())

    yield title, "".join(lines).replace(" ", "").replace("\r", "").upper()

REGEX_GZIPPED = re.compile(r'^.+\.gz$')
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
    print("Hi I am worker {} with PID {}".format(current_process(),getpid()))
    sys.stdin.flush()
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
    iter_keys = seqs.keys()
    kmer_freq = {}
    for seq_id in iter_keys:
        seq = seqs[seq_id].replace('-','')
        for idx, (kmername, kmer_seq, is_revcomp) in automaton.iter(seq):
            res.append((kmername, kmer_seq, is_revcomp, seq_id, idx))
    columns = ['kmername', 'seq', 'is_revcomp', 'contig_id', 'match_index']
    return pd.DataFrame(res, columns=columns)

def parallel_query_fasta_files(input_genomes,
                            automaton: Automaton,
                           n_threads: int = 1):
    results = []
    sys.stdin.flush()
    if n_threads == 1:

        for i in range(0,len(input_genomes)):
            stime = time.time()
            df = find_in_fasta(automaton,input_genomes[i])
            df['file_path'] = input_genomes[i]
            sample = os.path.basename(input_genomes[i])
            sample = re.sub(r"(\.fa$)|(\.fas$)|(\.fasta$)|(\.fna$)", "", sample)
            df['sample'] = sample
            df['is_pos_kmer'] = ~df.kmername.str.contains('negative')
            refpositions = [x for x, y in df.kmername.str.split('-')]
            df['refposition'] = [int(x.replace('negative', '')) for x in refpositions]
            results.append(df)

            print("{}\t{}".format(input_genomes[i],time.time()-stime))
    else:
        print("creating pool")
        pool = Pool(processes=n_threads)
        res = []
        for i in range(0, len(input_genomes)):
            res.append(pool.apply_async(find_in_fasta,  ( automaton, input_genomes[i]  )))
        print("closing pool")
        #cleanup
        pool.close()
        pool.join()

        for i in range(0,len(res)):
            df = res[i].get()
            df['file_path'] = input_genomes[i]
            sample = os.path.basename(input_genomes[i])
            sample = re.sub(r"(\.fa$)|(\.fas$)|(\.fasta$)|(\.fna$)", "", sample)
            df['sample'] = sample
            df['is_pos_kmer'] = ~df.kmername.str.contains('negative')
            refpositions = [x for x, y in df.kmername.str.split('-')]
            df['refposition'] = [int(x.replace('negative', '')) for x in refpositions]
            results.append(df)

    return pd.concat(results)

def parallel_query_contigs(input_genomes,
                            automaton: Automaton,
                           n_threads: int = 1):
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
    kmer_freq = {}
    for kmer_seq, freq in kmer_seq_counts.items():
        kmername, sequence, _ = automaton.get(kmer_seq)
        if not kmername in kmer_freq:
            kmer_freq[kmername] = 0
        kmer_freq[kmername]+= freq
        res.append((kmername, kmer_seq, freq))
    for kmername,freq in kmer_freq:
        res.append([kmername,freq])
    return pd.DataFrame(res, columns=['kmername','freq'])

def parallel_fastq_query(automaton: Automaton, fastqs, n_threads=1):
    results = []
    df = pd.DataFrame()
    if n_threads == 1:
        for file in fastqs:
            print(file)
            tmp = find_in_fastqs(automaton,file)
            tmp['file'] = file
            sample = os.path.basename(file)
            sample = re.sub(r"(\_1.fq$)|(\_2.fq$)|(\_R1.fq$)|(\_R2.fq$)|(\_1.fastq$)|(\_2.fastq$)|(\_R1.fastq$)|(\_R2.fastq$)|(\.fq$)|(\.fastq$)", "", sample)
            tmp['sample'] = sample
            results.append(tmp)
    else:

        pool = Pool(processes=n_threads)
        res = []
        for file in fastqs:
            res.append(pool.apply_async(find_in_fastqs,  ( automaton,file  )))
        #cleanup
        pool.close()
        pool.join()
        for i in range(0,len(res)):
            tmp = res[i].get()
            tmp['file'] = fastqs[i]
            sample = os.path.basename(fastqs[i])
            sample = re.sub(r"(\_1.fq$)|(\_2.fq$)|(\_R1.fq$)|(\_R2.fq$)|(\_1.fastq$)|(\_2.fastq$)|(\_R1.fastq$)|(\_R2.fastq$)|(\.fq$)|(\.fastq$)", "", sample)
            tmp['sample'] = sample
            results.append(tmp)


    df = pd.concat(results)
    refpositions = [x for x, y in df.kmername.str.split('-')]
    df['refposition'] = [int(x.replace('negative', '')) for x in refpositions]
    df['is_pos_kmer'] = ~df.kmername.str.contains('negative')

    return df