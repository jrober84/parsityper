#!/usr/bin/python
import logging
import sys
import time
from argparse import (ArgumentParser)
import os
import pandas as pd
from parsityper.version import __version__
from parsityper.helpers import init_console_logger
from multiprocessing import Pool
from parsityper.scheme import parseScheme, constructSchemeLookups


def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Parsityper scheme validator')
    parser.add_argument('--input_dir', type=str, required=True,
                        help='Typer result directory')
    parser.add_argument('--input_meta', type=str, required=True,
                        help='TSV file of sample_id, genotype')
    parser.add_argument('--scheme', type=str, required=False,
                        help='TSV formated kmer scheme',default='parsityper')
    parser.add_argument('--outdir', type=str, required=True,
                        help='output directory')
    parser.add_argument('--prefix', type=str, required=False,
                        help='output file prefix',default='parsityper')
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


def compare_profiles(profile_1,profile_2):
    num_features = len(profile_1)
    if num_features != len(profile_2):
        return -1
    dist = 0
    for i in range(0,num_features):
        f1 = profile_1[i]
        f2 = profile_2[i]
        if f1 == 0.5 or f2 == 0.5 or f1 == f2:
            continue
        dist+=1
    return dist

def pariwise_cmp_genotype(sample_id_1,profiles,i,num_genotypes,genotypes):
    issues = {}
    for k in range(i + 1, num_genotypes):
        sample_id_2 = genotypes[k]
        dist = compare_profiles(profiles[sample_id_1], profiles[sample_id_2])
        if dist == -1:
            logging.ERROR("Profiles for {} and {} do not have the same length".format(sample_id_1, sample_id_2))
        if dist == 0:
            if not sample_id_1 in issues:
                issues[sample_id_1] = []
            issues[sample_id_1].append(sample_id_2)
    return issues

def validate_genotype_rules(profiles,n_threads=1):
    num_genotypes = len(profiles)
    genotypes = sorted(list(profiles.keys()))
    issues = {}
    if n_threads > 1:
        pool = Pool(processes=n_threads)
    results = []
    for i in range(0,num_genotypes):
        sample_id_1 = genotypes[i]
        if n_threads == 1:
            issues.update(pariwise_cmp_genotype(sample_id_1, profiles, i, num_genotypes, genotypes))
        else:
            results.append(pool.apply_async(pariwise_cmp_genotype,(sample_id_1, profiles, i, num_genotypes, genotypes)))
    for i in range(0,len(results)):
        issues.update(results.get())
    return issues

def validate_typer_dir(files):
    result = 0
    for file in files:
        if not os.path.isfile(file):
            result = -1
            logging.ERROR("Required file {} is not found, check prefix and directory and try again".format(file))
            continue
        if os.path.getsize(file) == 0:
            result = -1
            logging.ERROR("Required file {} is empty, check prefix and directory and try again".format(file))
            continue
    return result

def calc_genotype_metrics(ground_truth,sample_assignments):
    metrics = {}
    total = 0
    for sample_id in ground_truth:
        if sample_id not in sample_assignments:
            continue
        genotype = ground_truth[sample_id]
        if not genotype in metrics:
            metrics[genotype] = {
                'total':0,
                'ambig':0,
                'conservative_TP':0,
                'conservative_FP': 0,
                'conservative_TN': 0,
                'conservative_FN': 0,
                'conservative_F1':0
            }
        metrics[genotype]['total'] += 1
        total+=1

    for sample_id in ground_truth:
        if sample_id not in sample_assignments:
            continue
        genotype = ground_truth[sample_id]

        classified_genotypes = sample_assignments[sample_id].strip().split(',')
        num_genotypes = len(classified_genotypes)
        if num_genotypes == 1 and genotype == classified_genotypes[0]:
            metrics[genotype]['conservative_TP'] += 1
            continue
        if num_genotypes == 1 and genotype != classified_genotypes[0]:
            metrics[genotype]['conservative_FN'] += 1
            continue

        for g in classified_genotypes:
            if g not in metrics:
                continue
            if g != genotype:
                metrics[g]['conservative_FP'] += 1

    for genotype in metrics:
        precision = metrics[genotype]['conservative_TP'] / (metrics[genotype]['conservative_TP'] + metrics[g]['conservative_FP'])
        recall = metrics[genotype]['conservative_TP'] / (metrics[genotype]['conservative_TP'] + metrics[genotype]['conservative_FN'])
        metrics[genotype]['conservative_F1'] = 2 * (precision - recall ) / (precision + recall )

    return metrics

def process_kmer_profile(df):
    samples = df.columns.tolist()
    uids = df.index.tolist()
    num_uids = len(uids)
    kmer_results = {}
    for sample_id in samples:





def run():
    cmd_args = parse_args()
    logger = init_console_logger(2)

    #input parameters
    prefix = cmd_args.prefix
    min_cov = cmd_args.min_cov
    min_cov_frac = cmd_args.min_cov_frac
    input_dir = cmd_args.input_dir
    input_meta = cmd_args.input_meta
    scheme_file = cmd_args.scheme
    prefix = cmd_args.prefix
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


    scheme = parseScheme(scheme_file)
    logger.info("Initializing scheme data structure from {}".format(scheme_file))
    scheme_info = constructSchemeLookups(scheme)

    #Identify genotypes which cannot be identified unambiguously from the scheme
    genotype_collision_report = os.path.join(outdir, "{}.genotype.collisions.txt".format(prefix))
    collisions = validate_genotype_rules(scheme_info['kmer_profiles'],n_threads)
    fh = open(genotype_collision_report,'w')
    out_str = []
    for genotype in collisions:
        row = [genotype] + collisions[genotype]
        out_str.append("\t".join([str(x) for x in row]))
    fh.write("\n".join(out_str))
    fh.close()
    del(collisions)
    sys.exit()
    #validate typer files
    profile_file = os.path.join(input_dir, "{}_analysis.sample.kmer.profiles.txt".format(prefix))
    typer_file  = os.path.join(input_dir,"{}_analysis.sample_composition.report.summary.txt".format(prefix))
    typer_files = [
        typer_file,
        profile_file,
        input_meta
    ]

    if validate_typer_dir(typer_files) == -1:
        logging.ERROR("Something went wrong with reading the typer files, please check the log messages")
        sys.exit()

    profile_df = pd.read_csv(profile_file,sep="\t",header=0,index_col=0)
    meta_df = pd.read_csv(input_meta,sep="\t",header=0)
    ground_truth = dict(zip(meta_df.sample_id, meta_df.genotype))
    typer_df = pd.read_csv(profile_file,sep="\t",header=0)
    sample_assignments = dict(zip(typer_df.sample_id, typer_df.compatible_genotypes))
    metrics = calc_genotype_metrics(ground_truth, sample_assignments)
    print(metrics)