#!/usr/bin/python
import time
from argparse import (ArgumentParser, FileType)
import logging, os, sys, re, collections, operator, math, shutil, datetime
from collections import Counter
import pandas as pd
from parsityper.helpers import validate_args, init_console_logger, read_tsv, scheme_to_biohansel_fasta,filter_biohansel_kmer_df, process_biohansel_kmer, \
    init_kmer_targets,generate_biohansel_kmer_names, get_scheme_template, generate_random_phrase, calc_md5, get_kmer_groups, get_kmer_group_mapping, get_sequence_files
import copy
from parsityper.bio_hansel import bio_hansel
import statistics
from parsityper.constants import PRIMER_SCHEMES, TYPING_SCHEMES
from parsityper.helpers import  profile_pairwise_distmatrix
from parsityper.visualizations import dendrogram_visualization


def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='SARS-COV-2 Kmer-analysis')
    parser.add_argument('--data_dir', type=str, required=True,
                        help='directory of sequence data for developing/validating typing rules')
    parser.add_argument('--input_meta', type=str, required=True,
                        help='TSV file of sample_id,genotype')
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
    parser.add_argument('--n_threads', type=str, required=False,
                        help='output directory',default=1)
    parser.add_argument('--list', type=str, required=False,
                        help='list in-built primer and typing schemes')
    return parser.parse_args()


def run():
    cmd_args = parse_args()
    logger = init_console_logger(2)
    is_args_ok = validate_args(cmd_args,logger)
    if not is_args_ok:
        logger.error("One or more command line arguments has an issue, please check the log messages and try again")
        sys.exit()
    return

    #input parameters
    prefix = cmd_args.prefix
    min_cov = cmd_args.min_cov
    min_cov_frac = cmd_args.min_cov_frac
