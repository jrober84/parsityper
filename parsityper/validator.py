#!/usr/bin/python
import time
from argparse import (ArgumentParser)
import os
import pandas as pd
from parsityper.version import __version__
from parsityper.helpers import init_console_logger, read_tsv, process_biohansel_kmer, \
    init_kmer_targets,get_kmer_groups, get_kmer_group_mapping, summarize_kmer_targets, read_samples, read_fasta
from parsityper.kmerSearch.kmerSearch import init_automaton_dict,perform_kmerSearch_fasta,perform_kmerSearch_fastq
from parsityper.helpers import  get_expanded_kmer_number
import multiprocessing as mp
from multiprocessing import pool
from parsityper.scheme import parseScheme, constructSchemeLookups

if mp.get_start_method(allow_none=True) != 'spawn':
        mp.set_start_method('spawn', force=True)

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='SARS-COV-2 Kmer-analysis')
    parser.add_argument('--input', type=str, required=True,
                        help='TSV file of sample_id,genotype,file_1,file_2 ')
    parser.add_argument('--scheme', type=str, required=False,
                        help='TSV formated kmer scheme',default='SARS-COV-2')
    parser.add_argument('--outdir', type=str, required=True,
                        help='output directory')
    parser.add_argument('--prefix', type=str, required=False,
                        help='output file prefix',default='SARS-COV-2_kmer_analysis')
    parser.add_argument('--min_cov', type=int, required=False,
                        help='Absolute minimum kmer coverage for fastq read detection default=auto determine coverage',default=50)
    parser.add_argument('--min_cov_frac', type=float, required=False,
                        help='Minimum percentage of total pool required for k-mer detection range 0 - 1.0 (default=0.05)', default=0.05)
    parser.add_argument('--max_frac_missing', type=float, required=False,
                        help='Maximum number of sequences allowed to be missing a target for it to stay included', default=0.25)
    parser.add_argument('--n_threads', type=int, required=False,
                        help='output directory',default=1)
    parser.add_argument('-V', '--version', action='version', version='%(prog)s {}'.format(__version__))
    return parser.parse_args()



def run():
    cmd_args = parse_args()
    logger = init_console_logger(2)

    #input parameters
    prefix = cmd_args.prefix
    min_cov = cmd_args.min_cov
    min_cov_frac = cmd_args.min_cov_frac
    input = cmd_args.input
    scheme_file = cmd_args.scheme
    prefix = cmd_args.prefix
    nthreads = cmd_args.n_threads
    outdir = cmd_args.outdir
    detection_limit = min_cov
    logger = init_console_logger(2)
    scheme_dict = {}
    scheme_df = read_tsv(scheme_file)
    genotypes = []

    # initialize analysis directory
    if not os.path.isdir(outdir):
        logger.info("Creating analysis results directory {}".format(outdir))
        os.mkdir(outdir, 0o755)
    else:
        logger.info("Results directory {} already exits, will overwrite any results files here".format(outdir))



    scheme = parseScheme(scheme_file)
    logger.info("Initializing scheme data structure from {}".format(scheme_file))
    scheme_info = constructSchemeLookups(scheme)

    # Init Ahocorasak automation objects
    logger.info("Initializing aho-corasick automation")
    aho = {'scheme': init_automaton_dict(scheme_info['uid_to_kseq'])}

    samples = read_samples(input)
    kmer_results = {}
    fasta_file_extensions = ['.fasta','.fas','.fa','.fna','.fasta.gz','.fas.gz','.fa.gz','.fna.gz']
    fastq_file_extensions = ['.fastq','.fq','.fastq.gz','.fq.gz']
    for sampleID in samples:
        read_set = []
        if samples[sampleID]['file_1'] is not None:
            read_set.append(samples[sampleID]['file_1'])
        if samples[sampleID]['file_2'] is not None:
            read_set.append(samples[sampleID]['file_2'])
        fileType = None
        for ext in fasta_file_extensions:
            if ext in read_set[0]:
                fileType = 'fasta'
                break
        for ext in fastq_file_extensions:
            if ext in read_set[0]:
                fileType = 'fastq'
                break

        if nthreads == 1:
            if fileType == 'fastq':
                kmer_results[sampleID] = perform_kmerSearch_fastq(scheme_info['uid_to_kseq'], scheme_info['kseq_to_uids'],
                                                                  aho['scheme'], read_set)
            else:
                kmer_results[sampleID] = perform_kmerSearch_fasta(scheme_info['uid_to_kseq'], scheme_info['kseq_to_uids'],
                                                                  aho['scheme'], read_fasta(read_set[0]), min_cov)
        else:
            if fileType == 'fastq':
                kmer_results[sampleID] = pool.apply_async(perform_kmerSearch_fastq, (
                scheme_info['uid_to_kseq'], scheme_info['kseq_to_uids'], aho['scheme'], read_set))
            else:
                kmer_results[sampleID] = pool.apply_async(perform_kmerSearch_fasta,
                                                          (scheme_info['uid_to_kseq'], scheme_info['kseq_to_uids'],
                                                           aho['scheme'], read_fasta(read_set[0]), min_cov))


