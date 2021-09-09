#!/usr/bin/python
import time
from argparse import (ArgumentParser, FileType)
import logging, os, sys, re, collections, operator, math, shutil, datetime, copy, statistics
from collections import Counter
import pandas as pd
from parsityper.helpers import process_rules, init_console_logger, read_tsv, get_valid_targets, process_biohansel_kmer, \
    init_kmer_targets,summarize_kmer_targets,get_kmer_groups, get_kmer_group_mapping, get_kmer_freq_by_genotype, evaluate_rules
from parsityper.kmerSearch import init_automaton_dict,parallel_query_fasta_files, parallel_fastq_query
from parsityper.typer import calc_kmer_ratio, identify_compatible_types, calc_type_coverage, type_occamization
from parsityper.validator import construct_genotype_profiles, identify_degenerate_kmers, read_samples
import multiprocessing as mp

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
    parser.add_argument('--min_partial_frac', type=float, required=False,
                        help='Minimum fraction of isolates positive for mutation for it to be partial 0 - 1.0 (default=0.1)', default=0.1)
    parser.add_argument('--min_positive_frac', type=float, required=False,
                        help='Minimum fraction of isolates positive for mutation for it to be positive 0 - 1.0 (default=0.1)', default=0.9)
    parser.add_argument('--min_positive_freq', type=int, required=False,
                        help='Minimum number of isolates positive for mutation for it to be valid for the scheme (default=1)', default=1)
    parser.add_argument('--max_frac_missing', type=float, required=False,
                        help='Maximum number of sequences allowed to be missing a target for it to stay included', default=0.25)
    parser.add_argument('--n_threads', type=int, required=False,
                        help='output directory',default=1)
    return parser.parse_args()

def run():
    print(mp.get_context())
    cmd_args = parse_args()
    logger = init_console_logger(2)

    #input parameters
    prefix = cmd_args.prefix
    min_cov = cmd_args.min_cov
    min_cov_frac = cmd_args.min_cov_frac
    input = cmd_args.input
    min_positive_frac = cmd_args.min_positive_frac
    min_partial_frac = cmd_args.min_partial_frac
    min_pos_freq = cmd_args.min_positive_freq

    scheme_file = cmd_args.scheme
    prefix = cmd_args.prefix
    max_frac_missing = cmd_args.max_frac_missing
    n_threads = cmd_args.n_threads
    outdir = cmd_args.outdir
    detection_limit = min_cov
    logger = init_console_logger(2)

    # initialize analysis directory
    if not os.path.isdir(outdir):
        logger.info("Creating analysis results directory {}".format(outdir))
        os.mkdir(outdir, 0o755)
    else:
        logger.info("Results directory {} already exits, will overwrite any results files here".format(outdir))

    scheme_dict = {}
    scheme_df = read_tsv(scheme_file)

    genotypes = []
    for row in scheme_df.itertuples():
        id = row.key
        pos_kmer = row.positive
        neg_kmer = row.negative
        scheme_dict["{}-{}".format(id,id)] = pos_kmer
        scheme_dict["negative{}-{}".format(id,id)] = neg_kmer
        positive_genotypes = str(row.positive_seqs)
        if positive_genotypes == 'nan':
            positive_genotypes = []
        else:
            positive_genotypes = positive_genotypes.split(',')
            genotypes.extend(positive_genotypes)

        partial_genotypes = str(row.partial_positive_seqs)
        if partial_genotypes == 'nan':
            partial_genotypes = []
        else:
            partial_genotypes = partial_genotypes.split(',')
            genotypes.extend(partial_genotypes)
    genotypes = list(set(genotypes))
    genotypes.sort()
    logging.info("Read {} kmer targets".format(len(scheme_df)))
    scheme_kmer_target_info = init_kmer_targets(scheme_df)


    scheme_kmer_target_keys = list(scheme_kmer_target_info.keys())
    scheme_kmer_target_keys.sort()
    scheme_kmer_target_keys = [str(i) for i in scheme_kmer_target_keys]
    scheme_kmer_groups = get_kmer_groups(scheme_kmer_target_info)
    scheme_target_to_group_mapping = get_kmer_group_mapping(scheme_kmer_target_info)
    logging.info("Initializing kmer search object")
    A = init_automaton_dict(scheme_dict)
    samples = read_samples(input)

    logging.info("Found {} fasta samples and {} fastq samples".format(len(samples['fasta']),len(samples['fastq'])))
    mode = 'fasta'
    reported_genotypes = {}
    if len(samples['fasta']) > 0:
        input_genomes = []
        for sample in samples['fasta']:
            input_genomes.append(sample['file_1'])
            reported_genotypes[sample['sample_id']] = sample['genotype']
        logging.info(
            "Begining kmer search on {} samples".format(len(samples['fasta'])))
        kmer_df = parallel_query_fasta_files(input_genomes,A,n_threads)
        kmer_df['freq'] = min_cov
    if len(samples['fastq']) >0:
        mode = 'fastq'
        input_genomes = []
        for sample in samples['fastq']:
            input_genomes.append(sample['file_1'])
            reported_genotypes[sample['sample_id']] = sample['genotype']
            if len(sample['file_2']) > 0:
                input_genomes.append(sample['file_2'])
            logging.info(
                "Begining kmer search on {} samples".format(len(samples['fastq'])))
        kmer_df = parallel_fastq_query(A,input_genomes,  n_threads)
    logging.info("Processing k-mer results")
    sample_kmer_results = process_biohansel_kmer(scheme_kmer_groups, scheme_target_to_group_mapping,
                                                 scheme_kmer_target_info, kmer_df, min_cov)
    sample_kmer_data = calc_kmer_ratio(sample_kmer_results, scheme_kmer_target_keys, min_cov)

    # Get list of strains which are compatible with the kmer information in the scheme
    logging.info("Identifying compatible genotypes")
    sample_kmer_data = identify_compatible_types(scheme_df, sample_kmer_data, min_cov_frac,
                                                 detection_limit=detection_limit)
    logging.info("Identifying genotype profiles")
    genotype_profiles = construct_genotype_profiles(scheme_kmer_target_info, genotypes)
    kmer_data = summarize_kmer_targets(scheme_kmer_target_keys,sample_kmer_data,min_cov)
    valid_targets = get_valid_targets(kmer_data, len(sample_kmer_data), min_pos_freq, max_frac_missing)
    genotype_kmer_data = get_kmer_freq_by_genotype(sample_kmer_data,reported_genotypes,min_cov,min_cov_frac)
    rules = {}
    for genotype in genotype_kmer_data:
        targets = genotype_kmer_data[genotype]
        for target in targets:
            if not target in rules:
                rules[target] = {
                    'positive': [],
                    'partials': []
                }

            count_pos = len(targets[target]['positive'])
            if count_pos == 0:
                targets[target]['positive'] = []
            count_neg = len(targets[target]['negative'])
            if count_neg == 0:
                targets[target]['negative'] = []
            count_missing = len(targets[target]['missing'])
            if count_missing == 0:
                targets[target]['missing'] = []


            total = len(set( targets[target]['positive'] + targets[target]['negative'] + targets[target]['missing'] ))

            if total == 0:
                continue

            if count_pos / total >= min_positive_frac:
                rules[target]['positive'].append(genotype)
            elif count_pos / total >= min_partial_frac:
                rules[target]['partials'].append(genotype)
    logging.info("Calculating scheme score")
    scores = process_rules(sample_kmer_data, valid_targets, rules, genotypes, reported_genotypes, min_cov, min_cov_frac)
    print(scores)
    scheme_tuning_iter1 = evaluate_rules(scores, rules, 0.05)

    logging.info("Iteration 1: {}".format(scheme_tuning_iter1['scheme_score']))
    logging.info("Calculating scheme score")

    scores = process_rules(sample_kmer_data, valid_targets, scheme_tuning_iter1['rules'], genotypes, reported_genotypes, min_cov, min_cov_frac)
    scheme_tuning_iter2 = evaluate_rules(scores, rules, 0.05)
    logging.info("Iteration 2: {}".format(scheme_tuning_iter2['scheme_score']))


    if scheme_tuning_iter1['scheme_score'] >= scheme_tuning_iter2['scheme_score']:
        scheme_rules = scheme_tuning_iter1['rules']
    else:
        scheme_rules = scheme_tuning_iter2['rules']

    for target in scheme_kmer_target_info:
        if not target in scheme_rules:
            scheme_kmer_target_info[target]['positive_seqs'] = ''
            scheme_kmer_target_info[target]['partial_positive_seqs'] = ''
        else:
            scheme_kmer_target_info[target]['positive_seqs'] = ', '.join(scheme_rules[target]['positive'])
            scheme_kmer_target_info[target]['partial_positive_seqs'] = ', '.join(scheme_rules[target]['partials'])
    pd.DataFrame().from_dict(scheme_kmer_target_info,orient='index').to_csv(os.path.join(outdir,"{}-updated.scheme.txt".format(prefix)),sep="\t",index=False)

if __name__ == '__main__':
    mp.freeze_support()
    run()