#!/usr/bin/python
import time
from argparse import (ArgumentParser, FileType)
import logging, os, sys, re, collections, operator, math, shutil, datetime
import pandas as pd
from parsityper.helpers import init_console_logger, read_tsv, process_biohansel_kmer, \
    init_kmer_targets,get_kmer_groups, get_kmer_group_mapping, summarize_kmer_targets
from parsityper.kmerSearch import init_automaton_dict,parallel_query_fasta_files, parallel_fastq_query
from parsityper.typer import calc_kmer_ratio, identify_compatible_types, calc_type_coverage, type_occamization
from parsityper.helpers import  get_expanded_kmer_number
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
    parser.add_argument('--max_frac_missing', type=float, required=False,
                        help='Maximum number of sequences allowed to be missing a target for it to stay included', default=0.25)
    parser.add_argument('--n_threads', type=int, required=False,
                        help='output directory',default=1)
    return parser.parse_args()

def identify_degenerate_kmers(sample_kmer_data,min_cov):
    kmer_data = {}
    for sample_id in sample_kmer_data:
        counts = sample_kmer_data[sample_id]['counts']
        for target in counts:
            if not target in kmer_data:
                kmer_data[target] = {
                    'unique': [],
                    'mixed': [],
                    'multi':[]
                }
            positive_count = counts[target]['positive']
            negative_count = counts[target]['negative']
            total = positive_count + negative_count
            if positive_count > min_cov or negative_count > min_cov :
                kmer_data[target]['multi'].append(sample_id)
            if positive_count > min_cov and negative_count > min_cov:
                kmer_data[target]['mixed'].append(sample_id)
            if total == min_cov:
                kmer_data[target]['unique'].append(sample_id)
    return kmer_data



def construct_genotype_profiles(scheme,genotypes):
    profiles = {}
    for genotype in genotypes:
        profiles[genotype] = {}

    informative_targets = []
    for target in scheme:
        positive_genotypes = str(scheme[target]['positive_seqs'])
        if positive_genotypes == 'nan':
            positive_genotypes = []
        else:
            positive_genotypes = positive_genotypes.split(',')
            informative_targets.append(int(target))

        partial_genotypes = str(scheme[target]['partial_positive_seqs'])
        if partial_genotypes == 'nan':
            partial_genotypes = []
        else:
            partial_genotypes = partial_genotypes.split(',')
            informative_targets.append(int(target))
        if not int(target) in informative_targets:
            continue
        for genotype in profiles:

            if genotype in positive_genotypes:
                profiles[genotype][target] = '+'
            elif genotype in partial_genotypes:
                profiles[genotype][target] = '?'
            else:
                profiles[genotype][target] = '-'
    return profiles

def identify_collisions(profiles,file):
    genotypes = list(profiles.keys())
    num_genotypes = len(genotypes)
    fh = open(file,'w')
    for i in range(0,num_genotypes):
        for k in range(i+1,num_genotypes):
            match = []
            mismatch = []
            ambig = []
            for target in profiles[genotypes[i]]:
                value_1 = profiles[genotypes[i]][target]
                value_2 = profiles[genotypes[k]][target]
                if value_1 == '-' and value_2 == '-':
                    match.append(target)
                elif value_1 == '+' and value_2 == '+':
                    match.append(target)
                elif value_1 == '-' and value_2 == '+':
                    mismatch.append(target)
                elif value_1 == '+' and value_2 == '-':
                    mismatch.append(target)
                elif value_1 == '?' and value_2 == '?':
                    match.append(target)
                elif value_1 == '?' and value_2 != '?':
                    ambig.append(target)
                elif value_1 != '?' and value_2 == '?':
                    ambig.append(target)
            type = 'ambig'
            if len(mismatch) > 0:
                type = 'distinct'
            elif len(ambig) == 0 and len(mismatch) ==0:
                type = 'indistinct'
            fh.write("{}\t{}\t{}\t{}\n".format(genotypes[i],genotypes[k],type,", ".join(mismatch)))
    fh.close()




def read_samples(file):
    df = pd.read_csv(file,header=0,sep="\t")
    samples = {
        'fasta':[],
        'fastq':[]
    }
    for row in df.itertuples():
        sample_id = row.sample_id
        genotype = row.genotype
        file_1 = row.file_1
        file_2 = row.file_2
        mode = 'fasta'
        if 'fastq' in file_1 or 'fq' in file_1:
            mode = 'fastq'
        samples[mode].append({'sample_id':sample_id,'genotype':genotype,'file_1':file_1,'file_2':file_2})
    return samples


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
    n_threads = cmd_args.n_threads
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
    scheme_kmer_target_info = init_kmer_targets(scheme_df)
    #identify_collisions(construct_genotype_profiles(scheme_kmer_target_info, genotypes),os.path.join(outdir,"{}-scheme-genotypes.txt".format(prefix)))
    scheme_kmer_target_keys = list(scheme_kmer_target_info.keys())
    count_exp_kmers = get_expanded_kmer_number(scheme_kmer_target_info)
    logger.info("Expected number of expanded kmers in this scheme is {}".format(count_exp_kmers))

    scheme_kmer_target_keys.sort()
    scheme_kmer_target_keys = [str(i) for i in scheme_kmer_target_keys]
    scheme_kmer_groups = get_kmer_groups(scheme_kmer_target_info)
    scheme_target_to_group_mapping = get_kmer_group_mapping(scheme_kmer_target_info)
    stime = time.time()
    A = init_automaton_dict(scheme_dict)
    samples = read_samples(input)

    mode = 'fasta'
    reported_genotypes = {}
    if len(samples['fasta']) > 0:
        input_genomes = []
        for sample in samples['fasta']:
            input_genomes.append(sample['file_1'])
            reported_genotypes[sample['sample_id']] = sample['genotype']
        kmer_df = parallel_query_fasta_files(input_genomes,A,n_threads)
        kmer_df['freq'] = min_cov
    if len(samples['fastq']) >0:
        mode = 'fastq'
        #stub need to fix
        input_genomes = []
        for sample in samples['fastq']:
            input_genomes.append(sample['file_1'])
            reported_genotypes[sample['sample_id']] = sample['genotype']
            if len(sample['file_2']) > 0:
                input_genomes.append(sample['file_2'])
        kmer_df = parallel_fastq_query(A,input_genomes,  n_threads)
        print(kmer_df)

    sample_kmer_results = process_biohansel_kmer(scheme_kmer_groups, scheme_target_to_group_mapping,
                                                 scheme_kmer_target_info, kmer_df, min_cov)
    sample_kmer_data = calc_kmer_ratio(sample_kmer_results, scheme_kmer_target_keys, min_cov)

    # Get list of strains which are compatible with the kmer information
    sample_kmer_data = identify_compatible_types(scheme_df, sample_kmer_data, min_cov_frac,
                                                 detection_limit=detection_limit)

    sample_kmer_data = calc_type_coverage(sample_kmer_data, scheme_df, min_cov_frac=min_cov_frac, min_cov=min_cov)
    sample_kmer_data = type_occamization(sample_kmer_data, scheme_df, min_cov_frac=min_cov_frac, min_cov=min_cov)
    sample_kmer_data = calc_type_coverage(sample_kmer_data, scheme_df, min_cov_frac=min_cov_frac, min_cov=min_cov,
                                          recalc=True)
    genotype_data = {}
    results = ["sample_id\treported\tgenotype,status\tpredicted_genotypes"]

    for sample in sample_kmer_data:
        reported_genotype = reported_genotypes[sample]
        if not reported_genotype in genotype_data:
            genotype_data[reported_genotype] = {'perfect':[],'ambig':[],'incorrect':[]}
        found_genotypes = sample_kmer_data[sample]['genotypes']['candidate_data']
        is_genotype_found = False
        if reported_genotype in found_genotypes:
            is_genotype_found = True
        if is_genotype_found and len(found_genotypes) == 1:
            status = 'Perfect Match'
            genotype_data[reported_genotype]['perfect'].append(sample)
        elif is_genotype_found and len(found_genotypes) > 1:
            status = 'Ambiguious Match'
            genotype_data[reported_genotype]['ambig'].append(sample)
        else:
            status = 'Incorrect Match'
            genotype_data[reported_genotype]['incorrect'].append(sample)
        out = [sample,reported_genotype,status,", ".join(found_genotypes)]
        results.append("\t".join(out))


    fh = open(os.path.join(outdir,"{}-sample.report.txt".format(prefix)),'w')
    fh.write("\n".join(results))
    fh.close()

    results = ["genotype\tcount\tperfect\tambig\tincorrect\tperc_perfect"]
    for genotype in genotype_data:
        count_perfect = len(genotype_data[genotype]['perfect'])
        count_ambig = len(genotype_data[genotype]['ambig'])
        count_incorrect = len(genotype_data[genotype]['incorrect'])
        total = count_perfect + count_ambig + count_incorrect
        results.append("\t".join(str(x) for x in [genotype,total,count_perfect,count_ambig,count_incorrect,count_perfect/total]))

    fh = open(os.path.join(outdir,"{}-genotype.report.txt".format(prefix)),'w')
    fh.write("\n".join(results))
    fh.close()


    kmer_target_data = summarize_kmer_targets(scheme_kmer_target_keys, sample_kmer_data, min_cov)

    results = ["kmer_id\tcount_found\tcount_missing\tcount_pos\tave_freq_pos\tcount_neg\tave_freq_neg\tmissing_samples"]
    for kmer_id in kmer_target_data:
        data = kmer_target_data[kmer_id]
        if len(data) == 0:
            data  = {
                'count_found':0,'count_positive':0,'count_negative':0,'missing_samples':list(sample_kmer_data.keys())
            }
        if not 'count_positive' in data:
            ave_freq_pos = 0
        elif len(data['count_positive']) > 0:
            ave_freq_pos = sum(data['count_positive']) / len(data['count_positive'])
        else:
            ave_freq_pos = 0
        if not 'count_negative' in data:
            ave_freq_neg = 0
        elif len(data['count_negative']) > 0:
            ave_freq_neg = sum(data['count_negative']) / len(data['count_negative'])
        else:
            ave_freq_neg = 0
        results.append("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(kmer_id,data['count_found'],
                                                           len(data['missing_samples']),
                                                           len(data['count_positive']),
                                                           ave_freq_pos,
                                                           len(data['count_negative']),
                                                           ave_freq_neg,
                                                           ",".join(data['missing_samples'])))
    fh = open(os.path.join(outdir,"{}-target.report.txt".format(prefix)),'w')
    fh.write("\n".join(results))
    fh.close()
    if mode == 'fasta':
        fh = open(os.path.join(outdir,"{}-target.report.txt".format(prefix)),'w')
        kmer_info = identify_degenerate_kmers(sample_kmer_data,min_cov)
        for target in kmer_info:
            out = [
                target,
                len(kmer_info[target]['unique']),
                len(kmer_info[target]['mixed']),
                len(kmer_info[target]['multi']),
                ", ".join(kmer_info[target]['multi'])
            ]
            fh.write("{}\n".format("\t".join([str(x) for x in out])))
        fh.close()

