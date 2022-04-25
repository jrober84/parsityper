#!/usr/bin/python
import datetime
import glob
import json
import logging
import os
import re
import statistics
import sys, time
from argparse import (ArgumentParser)
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
import plotly.express as px
from jinja2 import Template
from plotly.offline import plot

from parsityper.kmerSearch.kmerSearch import ksearch_fastq_lite,ksearch_fasta_lite,init_automaton_dict
from parsityper.helpers import read_fasta
from multiprocessing import Pool
from statistics import mean
from parsityper.constants import PRIMER_SCHEMES, TYPING_SCHEMES, TYPER_SAMPLE_SUMMARY_HEADER_BASE, FIGURE_CAPTIONS, \
    BATCH_HTML_REPORT
import re, glob, os
from parsityper.ext_tools.fastp import run_fastp
from parsityper.ext_tools.kmc import kmc_summary
from parsityper.ext_tools.lighter import run_lighter
from parsityper.ext_tools.jellyfish import run_jellyfish_count, parse_jellyfish_counts
from parsityper.helpers import read_fasta
from parsityper.helpers import validate_args, init_console_logger, read_tsv, \
    generate_random_phrase, calc_md5, nan_compatible_kmer_pairwise_distmatrix
from parsityper.kmerSearch.kmerSearch import init_automaton_dict, perform_kmerSearch_fasta, perform_kmerSearch_fastq, \
    process_kmer_results
from parsityper.scheme import parseScheme, constructSchemeLookups, detectAmbigGenotypes
from parsityper.version import __version__
from parsityper.visualizations import dendrogram_visualization, plot_mds, create_heatmap
import scipy
from scipy.cluster.hierarchy import fcluster
from scipy.cluster.hierarchy import dendrogram, linkage

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Parsimony based sample Kmer-analysis')
    parser.add_argument('--data_dir', type=str, required=False,
                        help='directory of fasta/fastq files')
    parser.add_argument('--se', type=str, required=False,
                        help='single-end fastq read')
    parser.add_argument('--R1', type=str, required=False,
                        help='paired-end fwd fastq read')
    parser.add_argument('--R2', type=str, required=False,
                        help='paired-end rev fastq read')
    parser.add_argument('--scheme', type=str, required=True,
                        help='TSV formated kmer scheme', default='')
    parser.add_argument('--scheme_meta', type=str, required=False,
                        help='Metadata of genotypes in supplied scheme')
    parser.add_argument('--sample_meta', type=str, required=False,
                        help='Metadata of samples to be processed')
    parser.add_argument('--outdir', type=str, required=True,
                        help='output directory')
    parser.add_argument('--prefix', type=str, required=False,
                        help='output file prefix', default='parsityper_analysis')
    parser.add_argument('--primers', type=str, required=False,
                        help='TSV formated primer file', default=None)
    parser.add_argument('--type', type=str, required=True,
                        help='Treat data as single or multi source')
    parser.add_argument('--seqTech', type=str, required=True,
                        help='Explicit sequencing technology: Illumina or Nanopore')
    parser.add_argument('--correct', required=False,
                        help='Perform read correction',action='store_true')
    parser.add_argument('--trim', required=False,
                        help='Perform read trimming',action='store_true')
    parser.add_argument('--dedup', required=False,
                        help='Perform read deduplication',action='store_true')
    parser.add_argument('--min_read_len', type=int, required=False,
                        help='Filter reads shorter than this limit',
                        default=0)
    parser.add_argument('--max_read_len', type=int, required=False,
                        help='Filter reads longer than this limit',
                        default=0)
    parser.add_argument('--min_genome_size', type=int, required=False,
                        help='Flag genomes which are too small',
                        default=0)
    parser.add_argument('--max_genome_size', type=int, required=False,
                        help='Flag genomes which are too large',
                        default=-1)
    parser.add_argument('--min_genome_cov_depth', type=int, required=False,
                        help='Flag samples where the raw sequencing depth is too low',
                        default=200)
    parser.add_argument('--min_cov', type=int, required=False,
                        help='Absolute minimum kmer coverage for fastq read detection default=auto determine coverage',
                        default=50)
    parser.add_argument('--min_cov_frac', type=float, required=False,
                        help='Minimum percentage of total pool required for k-mer detection range 0 - 1.0 (default=0.05)',
                        default=0.05)
    parser.add_argument('--max_mixed_sites', type=int, required=False,
                        help='Maximum number of sites allowed to have both kmer states', default=0)
    parser.add_argument('--max_missing_sites', type=int, required=False,
                        help='Maximum number of sites allowed to be missing', default=500)
    parser.add_argument('--genotype_dist_cutoff', type=float, required=False,
                        help='Maximum profile distance for considering a genotype result valid', default=0.0)
    parser.add_argument('--no_template_control', type=str, required=False,
                        help='SampleID of the negative control sample')
    parser.add_argument('--positive_control', type=str, required=False,
                        help='SampleID of the positive control sample')
    parser.add_argument('--no_plots', required=False,
                        help='suppress making plots, required for large datasets',action='store_true')
    parser.add_argument('--force_plots', required=False,
                        help='Try plotting for datasets > 1000',action='store_true')
    parser.add_argument('--delete_processed_reads', required=False,
                        help='Delete processed reads after completion',action='store_true')
    parser.add_argument('--typer_only', required=False,
                        help='Skip read preprocessing and stat calcuations', action='store_true')
    parser.add_argument('--max_features', type=str, required=False,
                        help='max gene features to report before collapsing', default=15)
    parser.add_argument('--n_threads', type=str, required=False,
                        help='output directory', default=1)
    parser.add_argument('-V', '--version', action='version', version='%(prog)s {}'.format(__version__))

    return parser.parse_args()

def checkPrimerKmerOverlap(scheme_info,primer_kmers):
    overlapping_kmers = {}
    for kSeq in scheme_info['kseq_to_uids']:
        for pName in primer_kmers:
            pSeq = primer_kmers[pName]
            if kSeq in pSeq:
                if not kSeq in overlapping_kmers:
                    overlapping_kmers[kSeq] = []
                overlapping_kmers[kSeq].append(pSeq)
    return overlapping_kmers

def perform_kmer_searching(seqManifest,scheme_info,min_cov,n_threads=1):
    aho = init_automaton_dict(scheme_info['uid_to_kseq'])
    num_kmers = len(scheme_info['uid_to_kseq'])
    kmer_counts = {}

    if n_threads > 1:
        pool = Pool(processes=n_threads)

    res = []
    fasta_files = []
    fastq_files = []
    for sample_id in seqManifest:
        fileType = seqManifest[sample_id]['file_type']
        read_set = seqManifest[sample_id]['raw_seq_files']
        if fileType == 'fasta':
            fasta_files.append(read_set[0])
        elif fileType == 'fastq':
            fastq_files.append(read_set)

    if len(fasta_files) > 0:
        res = []
        n = int(len(fasta_files) / n_threads)
        logging.info("Dividing fasta files into {} chunks of size {}".format(int(len(fasta_files) / n),n))
        batches = [fasta_files[i:i + n] for i in range(0, len(fasta_files), n)]
        for batch in batches:
            seqs = {}
            for file in batch:
                seqs.update(read_fasta(file))
            if n_threads > 1:
                res.append(pool.apply_async(ksearch_fasta_lite, (aho, num_kmers, seqs)))
            else:
                kmer_counts.update( ksearch_fasta_lite(aho, num_kmers, seqs))

    for sample_id in seqManifest:
        fileType = seqManifest[sample_id]['file_type']
        if len(seqManifest[sample_id]['processed_reads']) == 0:
            read_set = seqManifest[sample_id]['raw_seq_files']
        else:
            read_set = seqManifest[sample_id]['processed_reads']
        if fileType == 'fasta':
            continue
        if n_threads > 1:
            res.append( pool.apply_async(ksearch_fastq_lite, (aho, num_kmers, sample_id, read_set)))
        else:
            kmer_counts.update( ksearch_fastq_lite(aho, num_kmers, sample_id, read_set))

    if n_threads > 1:
        pool.close()
        pool.join()
        for result in res:
            result = result.get()
            for sample_id in result:
                kmer_counts[sample_id] = result[sample_id]

    for kseq in scheme_info['kseq_to_uids']:
        uids = scheme_info['kseq_to_uids'][kseq]
        for sample_id in kmer_counts:
            value = 0
            for uid in uids:
                value += kmer_counts[sample_id][uid]
            fileType = seqManifest[sample_id]['file_type']
            if fileType == 'fasta':
                value = value * min_cov
            for uid in uids:
                kmer_counts[sample_id][uid] = value

    return kmer_counts

def filter_low_freq_kmers(kmer_counts,num_kmers,min_freq=1):
    positions = range(0,num_kmers)
    for sample_id in kmer_counts:
        for i in positions:
            if kmer_counts[sample_id][i] < min_freq:
                kmer_counts[sample_id][i] = 0
    return kmer_counts

def filter_low_frac_kmers(kmer_counts,scheme_info,min_kmer_frac):
    sites = scheme_info['mutation_to_uid']
    for sample_id in kmer_counts:
        for site in sites:
            site_total_freq = 0
            for uid in sites[site]:
                site_total_freq+= kmer_counts[sample_id][uid]
            if site_total_freq == 0:
                continue
            for uid in sites[site]:
                kmer_frac = kmer_counts[sample_id][uid] / site_total_freq
                if kmer_frac < min_kmer_frac:
                    kmer_counts[sample_id][uid] = 0
    return kmer_counts


def filter_contam_kmers(kmer_counts, negative_control_id,scheme_info, min_kmer_frac):
    sites = scheme_info['mutation_to_uid']

    #init negative control
    neg_control = kmer_counts[negative_control_id]
    neg_control_sites = {}
    for site in sites:
        site_total_freq = 0
        for uid in sites[site]:
            site_total_freq += neg_control[uid]
        if site_total_freq == 0:
            continue
        neg_control_sites[site] = {'ref':[],'alt':[]}
        for uid in sites[site]:
            kmer_frac = neg_control[uid] / site_total_freq
            if kmer_frac < min_kmer_frac:
                neg_control[uid] = 0
            else:
                state = scheme_info['uid_to_state'][uid]
                neg_control_sites[site][state].append(uid)

    #subtract negative from each sample
    for sample_id in kmer_counts:
        for site in neg_control_sites:
            site_total_freq = 0

            for uid in sites[site]:
                site_total_freq += kmer_counts[sample_id][uid]

            if site_total_freq == 0:
                continue

            count_ref = 0
            count_alt = 0

            for uid in sites[site]:
                kmer_frac = kmer_counts[sample_id][uid] / site_total_freq
                if kmer_frac < min_kmer_frac:
                    kmer_counts[sample_id][uid] = 0
                else:
                    state = scheme_info['uid_to_state'][uid]
                    if state == 'ref':
                        count_ref+=1
                    else:
                        count_alt+=1

            if count_ref == 0 or count_alt == 0:
                continue

            for uid in sites[site]:
                kmer_counts[sample_id][uid] = 0

    return kmer_counts

def calc_average_kmer_freq(kmer_counts):
    averages = {}
    for sample_id in kmer_counts:
        averages[sample_id] = mean(kmer_counts[sample_id])
    return averages


def calc_site_frac(kmer_counts, scheme_info):
    sites = list(scheme_info['mutation_to_uid'].keys())
    num_sites = len(sites)
    site_ranges = range(0,num_sites)
    fracs = {}
    for sample_id in kmer_counts:
        fracs[sample_id] = [-1] * num_sites
        for i in site_ranges:
            site_total_freq = 0
            count_ref = 0
            count_alt = 0
            site = sites[i]
            for uid in scheme_info['mutation_to_uid'][site]:
                state = scheme_info['uid_to_state'][uid]
                site_total_freq += kmer_counts[sample_id][uid]
                if state == 'ref':
                    count_ref += kmer_counts[sample_id][uid]
                else:
                    count_alt += kmer_counts[sample_id][uid]
            if site_total_freq == 0:
                continue

            frac = count_alt / site_total_freq
            fracs[sample_id][i] = frac

    return fracs

def get_mixed_sites(kmer_fracs,scheme_info,min_frac):
    sites = list(scheme_info['mutation_to_uid'].keys())
    mixed_sites = {}
    for sample_id in kmer_fracs:
        mixed_sites[sample_id] = []
        for i in range(0,len(kmer_fracs[sample_id])):
            frac = kmer_fracs[sample_id][i]
            if frac >= min_frac and frac < 1:
                site = sites[i]
                mixed_sites[sample_id].append(site)
    return mixed_sites



def comp_mutation_profiles(mutation_fracs,scheme_info):
    sites = list(scheme_info['mutation_to_uid'].keys())
    num_mutations = len(sites)
    mut_range = range(0,num_mutations)
    genotype_profiles = scheme_info['mutation_profiles']
    dists = {}
    for sample_id in mutation_fracs:
        dists[sample_id] = {}
        for genotype in genotype_profiles:
            dist = 0
            total = 0
            for i in mut_range:
                rBase = mutation_fracs[sample_id][i]
                qBase = genotype_profiles[genotype][i]
                if qBase != 0.5 and rBase != -1:
                    total+=1
                if rBase == -1 or qBase == -1 or qBase == 0.5 or rBase == qBase or (rBase > 0 and rBase < 1):
                    continue
                dist+=1

            if total > 0:
                dists[sample_id][genotype] = dist / total
            else:
                dists[sample_id][genotype] = -1
            #print("{}\t{}\t{}".format(sample_id,genotype,dists[sample_id][genotype] ))
    return dists

def comp_kmer_profiles(kmer_counts,scheme_info,n_threads=1):
    kmer_list = list(scheme_info['uid_to_kseq'].keys())
    kmer_set = set(kmer_list)
    num_kmers = len(kmer_list)
    alt_kmers = scheme_info['inf_alt_uids']
    geno_rules = scheme_info['genotype_rule_sets']
    for genotype in geno_rules:
        for key in geno_rules[genotype]:
            geno_rules[genotype][key] = set(geno_rules[genotype][key])

    mutation_to_uid = scheme_info['mutation_to_uid']
    for mutation_key in mutation_to_uid:
        mutation_to_uid[mutation_key] = set(mutation_to_uid[mutation_key])

    genotype_inf_kmers = {}
    for genotype in geno_rules:
        genotype_inf_kmers[genotype] = geno_rules[genotype]['positive_uids'] | geno_rules[genotype]['partial_uids']

    if n_threads > 1:
        pool = Pool(processes=n_threads)
    res = []
    dists = {}
    for sample_id in kmer_counts:
        dists[sample_id] = {}
        exclude_sites = []
        detected_kmers = [i for i, v in enumerate( kmer_counts[sample_id]) if v > 0]
        detected_kmers = set(detected_kmers)

        for mutation_key in scheme_info['mutation_to_uid']:
            detected_mut_kmers = mutation_to_uid[mutation_key] & detected_kmers
            uids = mutation_to_uid[mutation_key]
            if len(detected_mut_kmers) == 0:
                exclude_sites += uids

        exclude_sites = set(exclude_sites)
        valid_kmers = kmer_set - exclude_sites
        if n_threads == 1:
            dists.update(calc_geno_dist(sample_id, detected_kmers, alt_kmers, valid_kmers, geno_rules, genotype_inf_kmers))
        else:
            res.append(pool.apply_async(calc_geno_dist,(sample_id, detected_kmers, alt_kmers, valid_kmers, geno_rules, genotype_inf_kmers)))

    if n_threads > 1:
        pool.close()
        pool.join()
        for result in res:
            result = result.get()
            for sample_id in result:
                dists[sample_id] = result[sample_id]

    return dists

def calc_geno_dist(sample_id,detected_kmers,alt_kmers,valid_kmers,geno_rules,genotype_inf_kmers):
    dists = {}
    dists[sample_id] ={}
    for genotype in geno_rules:
        dist = 0
        informative_uids = genotype_inf_kmers[genotype] & valid_kmers
        num_inf = len(informative_uids)
        matched = []
        if num_inf > 0:
            genotype_req_uids = geno_rules[genotype]['positive_uids'] & valid_kmers
            matched = detected_kmers & genotype_req_uids
            atypical = ((detected_kmers - informative_uids ) & alt_kmers)
            mismatched = (genotype_req_uids - matched) | atypical
            dist = len(mismatched) / (num_inf + len(atypical))
        if len(matched) == 0:
            dist = 1
        #print("{}\t{}\t{}".format(sample_id,genotype,dist))
        dists[sample_id][genotype] = dist

    return dists

def create_seq_manifest(input_dir):
    file_dict = find_seq_files(input_dir)
    seq_manifest = {}
    for sampleName in file_dict:
        seq_manifest[sampleName] = {
            'seq_type':'',
            'data_type':'sample',
            'is_valid':True,
            'is_compressed':False,
            'is_paired':False,
            'seq_files':[],
            'messages':[]
        }
        if len(file_dict[sampleName]) > 2:
            seq_manifest[sampleName]['messages'].append(
                "Warning {} files associated with the same sample id, only paired reads are supported for multiple files".formatlen(
                    file_dict[sampleName]))
            seq_manifest[sampleName]['is_valid'] = False
        if len(file_dict[sampleName]) == 2:
            seq_manifest[sampleName]['is_paired'] = True

        for i in range(0,len(file_dict[sampleName])):
            row = file_dict[sampleName][i]
            if i == 0:
                seq_manifest[sampleName]['seq_type'] = row['seq_type']
                seq_manifest[sampleName]['is_compressed'] = row['is_compressed']
                if row['contains_paired_info']:
                    seq_manifest[sampleName]['is_paired'] = True
            seq_manifest[sampleName]['seq_files'].append(row['file_path'])

    return seq_manifest

def find_seq_files(input_dir):
    file_dict = {}
    fasta_file_extensions = ['.fasta','.fas','.fa','.fna','.fasta.gz','.fas.gz','.fa.gz','.fna.gz']
    fastq_file_extensions = ['.fastq','.fq','.fastq.gz','.fq.gz']
    paired_end_notations = ['_R1','_R2','_1','_2',]
    if input_dir[-1] != '/':
        input_dir = "{}/".format(input_dir)
    file_list = glob.glob("{}**".format(input_dir), recursive=True)

    for file in file_list:
        fileType = ''
        fileName = os.path.basename(file)
        sampleName = fileName
        for ext in fasta_file_extensions:
            if re.search("{}$".format(ext),file):
                fileType = 'fasta'
                sampleName = re.sub(r"{}$".format(ext), '', sampleName )
                break

        if fileType == '':
            for ext in fastq_file_extensions:
                if re.search("{}$".format(ext),file):
                    fileType = 'fastq'
                    sampleName = re.sub(r"{}$".format(ext), '', sampleName )
                    break

        #skip files which do not contain these extections
        if fileType == '':
            continue

        is_compressed = False
        if file[-2:] == 'gz':
            is_compressed = True

        contains_paired_info = False


        for n in paired_end_notations:
            if sampleName != re.sub(r"{}$".format(n), '', sampleName):
                contains_paired_info = False
                sampleName = re.sub(r"{}$".format(n), '', sampleName)
                break

        entry = {
            'seq_type':fileType,
            'file_path':file,
            'sample_name':sampleName,
            'sequence_type':fileType,
            'is_compressed':is_compressed,
            'contains_paired_info':contains_paired_info,
            'seq_extension':ext,
            'paired_extension':''
        }

        if contains_paired_info:
            entry['paired_extension'] = n

        if not sampleName in file_dict:
            file_dict[sampleName] = []
        file_dict[sampleName].append(entry)

    return file_dict

def filter_dict_by_dist(samples,max_dist):
    filt = {}
    for sample_id in samples:
        if samples[sample_id] <= max_dist:
            filt[sample_id] = samples[sample_id]
    return filt


def assign_genotypes(kmer_counts,mutation_geno_dists,kmer_geno_dists,scheme_info,max_dist):
    kmer_list = list(scheme_info['uid_to_kseq'].keys())
    num_kmers = len(kmer_list)
    kmer_range = range(0, num_kmers)
    geno_rules = scheme_info['genotype_rule_sets']
    entropies = scheme_info['uid_to_entropy']
    dists = {}
    sample_genotypes = {}
    for sample_id in kmer_counts:
        dists[sample_id] = {}
        detected_kmers = []
        exclude_sites = []
        for i in kmer_range:
            if kmer_counts[sample_id][i] > 0:
                detected_kmers.append(i)

        detected_kmers = set(detected_kmers)
        for mutation_key in scheme_info['mutation_to_uid']:
            detected_mut_kmers = set(scheme_info['mutation_to_uid'][mutation_key]) & detected_kmers
            uids = scheme_info['mutation_to_uid'][mutation_key]
            if len(detected_mut_kmers) == 0:
                exclude_sites += uids
                continue

        exclude_sites = set(exclude_sites)
        valid_kmers = set(kmer_list) - exclude_sites

        mut_valid_genotypes = filter_dict_by_dist(mutation_geno_dists[sample_id], max_dist)
        kmer_valid_genotypes = filter_dict_by_dist(kmer_geno_dists[sample_id], max_dist)

        if len(kmer_valid_genotypes) == 0 and len(mut_valid_genotypes) == 0:
            sample_genotypes[sample_id] = {}
            continue
        elif  len(kmer_valid_genotypes) > 0 :
            candidate_genotypes = list(kmer_valid_genotypes.keys())
        else:
            candidate_genotypes = list(mut_valid_genotypes.keys())

        valid_kmer_entropies = {}
        for uid in valid_kmers:
            valid_kmer_entropies[uid] = entropies[uid]

        valid_kmer_entropies = {k: v for k, v in sorted(valid_kmer_entropies.items(), key=lambda item: item[1])}

        geno_kmer_assignments = {}
        geno_kmer_counts = {}
        genotyping_uids = []
        geno_ovl_kmer_count = {}
        for genotype in candidate_genotypes:
            geno_kmer_assignments[genotype] = []
            geno_kmer_counts[genotype] = 0
            genotyping_uids.extend(geno_rules[genotype]['positive_uids'])
            genotyping_uids.extend(geno_rules[genotype]['partial_uids'])
            geno_ovl_kmer_count[genotype] = len(detected_kmers & set(geno_rules[genotype]['positive_uids']))
        genotyping_uids = set(genotyping_uids) & set(list(valid_kmer_entropies.keys()))
        unassigned_kmers = genotyping_uids
        geno_ovl_kmer_count = {k: v for k, v in sorted(geno_ovl_kmer_count.items(), key=lambda item: item[1],reverse=True)}

        if len(candidate_genotypes) == 1:
            sample_genotypes[sample_id] = {genotype:list(genotyping_uids)}
            continue

        #use more specific kmer dists before mutations if needed
        dist_to_geno = {}
        for genotype in kmer_valid_genotypes:
            dist = str(kmer_valid_genotypes[genotype])
            if not dist in dist_to_geno:
                dist_to_geno[dist] = []
            dist_to_geno[dist].append(genotype)

        if len(dist_to_geno) == 0:
            for genotype in mut_valid_genotypes:
                dist = str(mut_valid_genotypes[genotype])
                if not dist in dist_to_geno:
                    dist_to_geno[dist] = []
                dist_to_geno[dist].append(genotype)

        #Assign singleton alt positive kmers first and then those that have different ditances
        assigned_kmers = []
        for uid in unassigned_kmers:
            state = scheme_info['uid_to_state'][uid]
            if state == 'ref':
                continue
            potential_genotypes = []
            for genotype in candidate_genotypes:
                if uid in geno_rules[genotype]['positive_alt']:
                    potential_genotypes.append(genotype)

            if len(potential_genotypes) == 1:
                geno_kmer_assignments[genotype].append(uid)
                geno_kmer_counts[genotype] += 1
                assigned_kmers.append(uid)
                continue

        unassigned_kmers = unassigned_kmers - set(assigned_kmers)
        geno_kmer_counts = {k: v for k, v in sorted(geno_kmer_counts.items(), key=lambda item: item[1],reverse=True)}

        # Assign singleton ref positive kmers
        for uid in unassigned_kmers:
            state = scheme_info['uid_to_state'][uid]
            if state == 'alt':
                continue
            potential_genotypes = []
            for genotype in geno_kmer_counts:
                if uid in geno_rules[genotype]['positive_ref']:
                    potential_genotypes.append(genotype)
            if len(potential_genotypes) == 1:
                geno_kmer_assignments[genotype].append(uid)
                geno_kmer_counts[genotype] += 1
                assigned_kmers.append(uid)
                continue

        unassigned_kmers = unassigned_kmers - set(assigned_kmers)

        geno_kmer_counts = {k: v for k, v in sorted(geno_kmer_counts.items(), key=lambda item: item[1], reverse=True)}

        grouped_counts = {}
        for genotype in geno_kmer_counts:
            count = geno_kmer_counts[genotype]
            if not count in grouped_counts:
                grouped_counts[count] = []
            grouped_counts[count].append(genotype)


        for count in grouped_counts:
            genotypes = grouped_counts[count]
            subset = {}
            for genotype in genotypes:
                subset[genotype] =  geno_ovl_kmer_count[genotype]
            subset = {k: v for k, v in sorted(subset.items(), key=lambda item: item[1], reverse=True)}

            for genotype in subset:
                ovl = unassigned_kmers & set(geno_rules[genotype]['positive_uids'])
                geno_kmer_counts[genotype] += len(ovl)
                geno_kmer_assignments[genotype].extend(list(ovl))
                assigned_kmers.extend(list(ovl))
                unassigned_kmers = unassigned_kmers - ovl

        geno_kmer_counts = {k: v for k, v in sorted(geno_kmer_counts.items(), key=lambda item: item[1], reverse=True)}


        # filter out genotypes without any assigned kmers
        filt = {}
        for genotype in geno_kmer_counts:
            if geno_kmer_counts[genotype] == 0:
                del (geno_kmer_assignments[genotype])
                continue

            filt[genotype] = geno_kmer_counts[genotype]
        geno_kmer_counts = filt

        unassigned_kmers = unassigned_kmers - set(assigned_kmers)

        for genotype in geno_kmer_counts:
            ovl = unassigned_kmers & set(geno_rules[genotype]['partial_uids'])
            geno_kmer_counts[genotype] += len(ovl)
            geno_kmer_assignments[genotype].extend(list(ovl))
            assigned_kmers.extend(list(ovl))
            unassigned_kmers = unassigned_kmers - ovl
        sample_genotypes[sample_id] = geno_kmer_assignments

    return sample_genotypes

def init_sample_manifest(seqManifest,scheme_info,scheme_name,analysis_date,sample_type,seqTech):
    sampleManifest = {}
    for sampleID in seqManifest:
        logging.info("Processing sample: {}".format(sampleID))
        sampleManifest[sampleID] = {
            'sample_id': sampleID,
            'scheme':scheme_name,
            'analysis_date':analysis_date,
            'reported_sample_type':sample_type,
            'sequencing_technology': seqTech,
            'file_type': seqManifest[sampleID]['seq_type'],
            'data_type':seqManifest[sampleID]['data_type'],
            'est_genome_size':0,
            'num_unique_kmers':0,
            'num_counted_unique_kmers':0,
            'num_seq_files': len(seqManifest[sampleID]['seq_files']),
            'raw_seq_files': seqManifest[sampleID]['seq_files'],
            'total_reads_pre':0,
            'total_reads_post': 0,
            'total_bases_pre': 0,
            'total_bases_post': 0,
            'read_mean_len_pre':0,
            'read_mean_len_post': 0,
            'read_gc_pre': 0,
            'read_gc_post': 0,
            'read_insert_size_peak': 0,
            'read_duplication_rate': 0,
            'estimated_genome_cov':0,
            'processed_reads':[],
            'detected_sample_type': 'single',
            'qc_messages':[],
            'compatible_genotypes':[],
            'total_scheme_kmers': len(scheme_info['uid_to_kseq']),
            'detected_scheme_kmers':[],
            'num_detected_scheme_kmers':0,
            'ave_scheme_kmers_freq': 0,
            'total_scheme_mutations': len(scheme_info['mutation_to_uid']),
            'detected_scheme_mutations': [],
            'num_detected_scheme_mutations':0,
            'detected_scheme_mixed_mutations': [],
            'num_detected_scheme_mixed_mutations': 0,
            'mutation_profile':[],
            'mutation_profile_md5':'',
            'mutation_profile_phrase':'',
            'primary_genotype':'',
            'primary_genotype_frac':'',
        }
    logging.info("Read {} samples into the manifest".format(len(sampleManifest)))
    return sampleManifest

def QA_results(sampleManifest,min_coverage_depth,max_missing_sites,min_genome_size=0,max_genome_size=-1,max_mixed_sites=0):
    for sample_id in sampleManifest:
        if sampleManifest[sample_id]['num_detected_scheme_mixed_mutations'] > max_mixed_sites:
            detected_sample_type = 'multi'
        else:
            detected_sample_type = 'single'
        sampleManifest[sample_id]['detected_sample_type'] = detected_sample_type
        reported_sample_type = sampleManifest[sample_id]['reported_sample_type']
        if detected_sample_type != reported_sample_type:
            sampleManifest[sample_id]['qc_messages'].append('Warning: sample type mismatch, reported:{} predicted:{}'.format(reported_sample_type,detected_sample_type))

        cov = sampleManifest[sample_id]['estimated_genome_cov']

        if cov < min_coverage_depth and sampleManifest[sample_id]['file_type'] == 'fastq':
            sampleManifest[sample_id]['qc_messages'].append('Fail: low sequencing coverage {}'.format(cov))
        else:
            sampleManifest[sample_id]['estimated_genome_cov'] = 1

        missing = sampleManifest[sample_id]['total_scheme_mutations'] - sampleManifest[sample_id]['num_detected_scheme_mutations']

        if missing > max_missing_sites:
            sampleManifest[sample_id]['qc_messages'].append(
                'Fail: high number of missing scheme targets {}'.format(missing))

        genomeSize = sampleManifest[sample_id]['est_genome_size']
        num_targets_found = sampleManifest[sample_id]['est_genome_size']
        if genomeSize < min_genome_size and genomeSize != 0:
            sampleManifest[sample_id]['qc_messages'].append(
                'Fail: genome size {} too small'.format(num_targets_found))
        elif genomeSize > max_genome_size and max_genome_size >= min_genome_size:
            sampleManifest[sample_id]['qc_messages'].append(
                'Fail: genome size {} too large'.format(num_targets_found))

    return sampleManifest

def process_reads(sampleManifest,perform_read_correction,read_dir,seqTech,fastp_dir,min_read_len,trim_front_bp,trim_tail_bp,perform_read_dedup,n_threads=1):
    for sampleID in sampleManifest:

        genomeSize = sampleManifest[sampleID]['est_genome_size']
        fileType = sampleManifest[sampleID]['file_type']
        seqFiles = sampleManifest[sampleID]['raw_seq_files']
        if fileType == 'fasta':
            continue

        # Perform read correction if requested
        if perform_read_correction:
            # initialize read directory
            sample_read_dir = os.path.join(read_dir, "{}".format(sampleID))
            if seqTech == 'illumina':
                run_lighter(seqFiles, sample_read_dir, 17, genomeSize, n_threads=n_threads)
                sampleManifest[sampleID]['processed_reads'] = glob.glob("{}/*.cor.*".format(sample_read_dir))

        # Perform read preprocessing using fastp
        if len(sampleManifest[sampleID]['processed_reads']) > 0:
            read_set = sampleManifest[sampleID]['processed_reads']
        else:
            read_set = sampleManifest[sampleID]['raw_seq_files']

        merge = False
        if len(read_set) == 2:
            merge = True

        fastp_results = run_fastp(read_set, fastp_dir, sampleID, min_read_len=min_read_len, trim_front_bp=trim_front_bp,
                                  trim_tail_bp=trim_tail_bp, report_only=False,
                                  dedup=perform_read_dedup, merge_reads=merge, n_threads=n_threads)

        sampleManifest[sampleID]['total_reads_pre'] = fastp_results['summary']["before_filtering"]["total_reads"]
        sampleManifest[sampleID]['total_bases_pre'] = fastp_results['summary']["before_filtering"]["total_bases"]
        sampleManifest[sampleID]['gc_pre'] = fastp_results['summary']["before_filtering"]["gc_content"]
        sampleManifest[sampleID]['total_reads_post'] = fastp_results['summary']["after_filtering"]["total_reads"]
        sampleManifest[sampleID]['total_bases_post'] = fastp_results['summary']["after_filtering"]["total_bases"]
        sampleManifest[sampleID]['gc_post'] = fastp_results['summary']["after_filtering"]["gc_content"]
        if "read1_mean_length" in fastp_results['summary']["before_filtering"]:
            sampleManifest[sampleID]['read_mean_len_pre'] = fastp_results['summary']["before_filtering"][
                "read1_mean_length"]
            sampleManifest[sampleID]['read_mean_len_post'] = fastp_results['summary']["after_filtering"][
                "read1_mean_length"]
        else:
            ave_len = fastp_results['summary']["before_filtering"]['total_bases'] / \
                      fastp_results['summary']["before_filtering"]["total_reads"]
            sampleManifest[sampleID]['read_mean_len_pre'] = ave_len
            ave_len = fastp_results['summary']["after_filtering"]['total_bases'] / \
                      fastp_results['summary']["after_filtering"]["total_reads"]
            sampleManifest[sampleID]['read_mean_len_post'] = ave_len
        if 'duplication' in fastp_results:
            sampleManifest[sampleID]['duplication_rate'] = fastp_results["duplication"]

        if genomeSize > 0:
            sampleManifest[sampleID]['estimated_genome_cov'] = sampleManifest[sampleID]['total_bases_pre'] / genomeSize


        if len(read_set) == 2:
            sampleManifest[sampleID]['read_mean_len_pre'] = (
            fastp_results['summary']["before_filtering"]["read1_mean_length"])
            sampleManifest[sampleID]['read_mean_len_post'] = (
            fastp_results['summary']["after_filtering"]["read1_mean_length"])
            pR1 = os.path.join(fastp_dir, "{}_1.fastq".format(sampleID))
            pR2 = os.path.join(fastp_dir, "{}_2.fastq".format(sampleID))
            pM = os.path.join(fastp_dir, "{}.merged.fastq".format(sampleID))
            sampleManifest[sampleID]['processed_reads'] = [pR1, pR2, pM]
    return sampleManifest

def calc_genome_size(sampleManifest,outdir,min_cov,kLen=21,n_threads=1):
    results = {}
    if n_threads > 1:
        pool = Pool(processes=n_threads)

    for sample_id in sampleManifest:
        cov = int(sampleManifest[sample_id]['ave_scheme_kmers_freq'])
        if cov < 1:
            cov = min_cov
        fileType = sampleManifest[sample_id]['file_type']
        if fileType == 'fasta':
            cov = 1
            seqs = read_fasta(sampleManifest[sample_id]['raw_seq_files'][0])
            if len(seqs) == 1:
                sampleManifest[sample_id]['est_genome_size'] = len(seqs[next(iter(seqs))].replace('-',''))
                sampleManifest[sample_id]['num_unique_kmers'] = ''
                sampleManifest[sample_id]['num_counted_unique_kmers'] = ''
                continue

        tmpFile = os.path.join(outdir, "{}.kmc.tmp1".format(sample_id))
        tmpDir = os.path.join(outdir, "{}.kmc.tmp2".format(sample_id))

        # get GenomeSize using KMC based on fwd read or fasta
        if n_threads == 1:
            results[sample_id] = kmc_summary(sampleManifest[sample_id]['raw_seq_files'][0], tmpFile,
                               tmpDir, fileType,
                               cov, kLen, n_threads)
        else:
            results[sample_id] = pool.apply_async(kmc_summary,
                                                           (sampleManifest[sample_id]['raw_seq_files'][0], tmpFile,
                               tmpDir, fileType,
                               cov, kLen, n_threads))

    if n_threads > 1:
        pool.close()
        pool.join()

        for sample_id in results:
            results[sample_id] = results[sample_id].get()

    for sample_id in results:
        for field in results[sample_id]:
            sampleManifest[sample_id][field] = results[sample_id][field]
        sampleManifest[sample_id]['estimated_genome_cov'] = sampleManifest[sample_id]['total_bases_pre'] / sampleManifest[sample_id]['est_genome_size']
    return sampleManifest

def add_genotype_info(sampleManifest,scheme_info,mutation_fracs,genotype_assignments):
    sites = list(scheme_info['mutation_to_uid'].keys())
    num_sites = len(sites)
    site_ranges = range(0,num_sites)

    for sample_id in genotype_assignments:
        genotypes = genotype_assignments[sample_id]
        sampleManifest[sample_id]['compatible_genotypes'] = list(genotypes.keys())
        if len(genotypes) == 0:
            continue
        elif len(genotypes) == 1:
            sampleManifest[sample_id]['primary_genotype'] = next(iter(genotypes))
            sampleManifest[sample_id]['primary_genotype_frac'] = 1.0
        else:
            genotype_fracs = {}
            for genotype in genotypes:
                mutations = []
                for uid in genotypes[genotype]:
                    mutations.append( scheme_info['uid_to_mutation'][uid])
                mutations = set(mutations)
                fracs = []
                for i in site_ranges:
                    site = sites[i]
                    if site in mutations:
                        frac = mutation_fracs[sample_id][i]
                        if frac == 0:
                            frac = 1
                        fracs.append(frac)
                ave_frac = 0
                if len(fracs) > 1:
                    ave_frac = statistics.mean(fracs)
                genotype_fracs[genotype] = ave_frac
            genotype_fracs = {k: v for k, v in sorted(genotype_fracs.items(), key=lambda item: item[1], reverse=True)}
            genotype = next(iter(genotype_fracs))
            sampleManifest[sample_id]['primary_genotype'] = genotype
            sampleManifest[sample_id]['primary_genotype_frac'] = genotype_fracs[genotype]

    return sampleManifest

def add_kmer_count_data(sampleManifest,mutation_fracs,scheme_info,kmer_counts):
    sites = list(scheme_info['mutation_to_uid'].keys())
    phrase_lookup = {}
    for sample_id in kmer_counts:
        detected_kmers_ids = []
        detected_kmers_vals = []
        for idx, val in enumerate(kmer_counts[sample_id]):
            if val != 0:
                detected_kmers_ids.append(idx)
                detected_kmers_vals.append(val)
        sampleManifest[sample_id]['detected_scheme_kmers'] = detected_kmers_ids
        sampleManifest[sample_id]['num_detected_scheme_kmers'] = len(detected_kmers_ids)
        if len(detected_kmers_vals) > 0:
            sampleManifest[sample_id]['ave_scheme_kmers_freq'] = statistics.mean(detected_kmers_vals)
        detected_mut = {}
        num_detected_mut = 0
        for idx,frac in enumerate(mutation_fracs[sample_id]):
            site = sites[idx]
            if frac == -1:
                continue
            num_detected_mut+=1
            detected_mut[site] = frac
        mutation_profile = []
        for site in detected_mut:

            if detected_mut[site] <= 0:
                continue
            elif detected_mut[site] > 0 and detected_mut[site] < 1:
                site = "{}*".format(site)
            mutation_profile.append(site)

        sampleManifest[sample_id]['mutation_profile'] = ",".join([str(x) for x in mutation_profile])
        sampleManifest[sample_id]['mutation_profile_md5'] = calc_md5(sampleManifest[sample_id]['mutation_profile'])
        if sampleManifest[sample_id]['mutation_profile_md5'] in phrase_lookup:
            sampleManifest[sample_id]['mutation_profile_phrase'] = phrase_lookup[sampleManifest[sample_id]['mutation_profile_md5']]
        else:
            phrase = ' '.join(generate_random_phrase())
            phrase_lookup[sampleManifest[sample_id]['mutation_profile_md5']] = phrase
            sampleManifest[sample_id]['mutation_profile_phrase'] = phrase
        sampleManifest[sample_id]['num_detected_scheme_mutations'] = num_detected_mut
        sampleManifest[sample_id]['detected_scheme_mutations'] = list(detected_mut.keys())
    return sampleManifest

def write_sample_summary_results(sampleManifest,scheme_info,out_file,sample_metadata,genotype_metadata,max_features=20):
    gene_features = scheme_info['gene_features']
    gene_features.sort()
    num_gene_features = len(gene_features)
    header = []
    if len(sample_metadata) > 0:
        key = next(sample_metadata)
        sample_metadata_fields = list(sample_metadata[key].keys())
        header += sample_metadata_fields
    header += TYPER_SAMPLE_SUMMARY_HEADER_BASE
    if len(genotype_metadata) > 0:
        key = next(genotype_metadata)
        geno_fields = list(genotype_metadata[key].keys())
        header += geno_fields

    if num_gene_features  <= max_features:
        for feature in gene_features:
            if feature != 'intergenic':
                header.append("dna_{}".format(feature))
                header.append("aa_{}".format(feature))
            else:
                header.append("dna_{}".format(feature))
    else:
        header.extend(['dna_cds','aa_cds','dna_intergenic'])
    fh = open(out_file,'w')
    fh.write("{}\n".format("\t".join(header)))
    #scheme lookups
    uid_to_dna_name = scheme_info['uid_to_dna_name']
    uid_to_gene_feature = scheme_info['uid_to_gene_feature']
    uid_to_state = scheme_info['uid_to_state']
    uid_to_aa_name = scheme_info['uid_to_aa_name']
    for sample_id in sampleManifest:
        sample_data = {}
        for field in header:
            if field in sampleManifest[sample_id]:
                sample_data[field] = sampleManifest[sample_id][field]
            else:
                sample_data[field] = ''
        kmer_uids = sampleManifest[sample_id]['detected_scheme_kmers']

        sample_features_dna = {}
        sample_features_aa = {}
        for feature in gene_features:
            sample_features_dna[feature] = []
            sample_data["dna_{}".format(feature)] = []
            if feature != 'intergenic':
                sample_data["aa_{}".format(feature)] = []
                sample_features_aa[feature] = []

        # only report alt state mutations
        # assign each mutation to a gene feature, collapse if there are too many features
        for uid in kmer_uids:
            state = uid_to_state[uid]
            if state == 'ref':
                continue
            feature = uid_to_gene_feature[uid]
            dna_name = uid_to_dna_name[uid]
            sample_features_dna[feature].append(dna_name)
            if feature != 'intergenic':
                aa_name = uid_to_aa_name[uid]
                sample_features_aa[feature].append(aa_name)

        if num_gene_features <= max_features:
            for feature in gene_features:
                sample_data["dna_{}".format(feature)] = sample_features_dna[feature]
                if feature != 'intergenic':
                    sample_data["aa_{}".format(feature)] = sample_features_aa[feature]
        else:
            cds_dna = []
            cds_aa = []
            for feature in gene_features:
                if feature != 'intergenic':
                    cds_dna.extend(sample_features_dna[feature])
                    cds_aa.extend(sample_features_aa[feature])
            sample_data['dna_cds'] = cds_dna
            sample_data['aa_cds'] = cds_aa
            sample_data['dna_intergenic'] = sample_features_dna['intergenic']
        primary_genotype = sample_data['primary_genotype']
        if primary_genotype in genotype_metadata:
            for field in genotype_metadata[primary_genotype]:
                sample_data[field] = genotype_metadata[primary_genotype][field]
        if sample_id in sample_metadata:
            for field in sample_metadata[sample_id]:
                sample_data[field] = sample_metadata[sample_id][field]
        row = []
        for field in header:
            value = ''
            if field in sample_data:
                value = sample_data[field]
            if isinstance(value,list):
                value = list(set(value))
                value.sort()
                value = ", ".join([str(x) for x in value])
            row.append(str(value))
        fh.write("{}\n".format("\t".join(row)))
    fh.close()

def create_run_summary(sampleManifest,kmer_counts,scheme_info):
    mutations = {}
    for mutation_key in scheme_info['mutation_to_uid']:
        mutations[mutation_key] =  {
            'sample_count':0,
            'mixed_sample_count':0,
            'samples': [],
            'mixed_samples': [],
            'mutation_freq':0,
            'kmer_freqs': {},
            'kmer_ave_freqs': {},
            'neg_control_freqs':{},
            'neg_control_ave_freqs': 0,
            'pos_control_freqs': {},
            'pos_control_ave_freqs': 0,
        }
        for uid in scheme_info['mutation_to_uid'][mutation_key]:
            mutations[mutation_key]['kmer_freqs'][uid] = []
            mutations[mutation_key]['kmer_ave_freqs'][uid] = 0
            mutations[mutation_key]['neg_control_freqs'][uid] = []
            mutations[mutation_key]['pos_control_freqs'][uid] = []

    for sample_id in sampleManifest:
        data_type = sampleManifest[sample_id]['data_type']
        mixed = sampleManifest[sample_id]['detected_scheme_mixed_mutations']
        for mutation_key in mixed:
            mutations[mutation_key]['mixed_samples'].append(sample_id)
        detected = sampleManifest[sample_id]['detected_scheme_mutations']
        for mutation_key in detected:
            mutations[mutation_key]['samples'].append(sample_id)
            uids = scheme_info['mutation_to_uid'][mutation_key]
            for uid in uids:
                freq = kmer_counts[sample_id][uid]
                if freq > 0:
                    mutations[mutation_key]['kmer_freqs'][uid].append(freq)
                    if data_type == 'negative':
                        mutations[mutation_key]['neg_control_freqs'][uid].append(freq)
                    elif data_type == 'positive':
                        mutations[mutation_key]['pos_control_freqs'][uid].append(freq)

    for mutation_key in mutations:
        mutations[mutation_key]['sample_count'] = len(mutations[mutation_key]['samples'])
        mutations[mutation_key]['mixed_sample_count'] = len(mutations[mutation_key]['mixed_samples'])
        for uid in mutations[mutation_key]['kmer_freqs']:
            if len(mutations[mutation_key]['kmer_freqs'][uid]) > 0:
                mutations[mutation_key]['kmer_ave_freqs'][uid] = statistics.mean(mutations[mutation_key]['kmer_freqs'][uid])
        freqs = []
        for uid in mutations[mutation_key]['kmer_ave_freqs']:
            freq = mutations[mutation_key]['kmer_ave_freqs'][uid]
            if freq > 0:
                freqs.append(freq)
        if len(freqs) > 0:
            mutations[mutation_key]['mutation_freq'] = statistics.mean(freqs)

        for uid in mutations[mutation_key]['kmer_freqs']:
            if len(mutations[mutation_key]['kmer_freqs'][uid]) > 0:
                mutations[mutation_key]['kmer_ave_freqs'][uid] = statistics.mean(mutations[mutation_key]['kmer_freqs'][uid])

    return mutations

def write_run_summary(run_summary,outdir,num_samples,prefix):
    summary_out = os.path.join(outdir,"{}.run.summary.txt")
    mutations_out = os.path.join(outdir,"{}.run.mutations.report.txt")
    kmers_out = os.path.join(outdir, "{}.run.mutations.report.txt")

    site_summary = {}
    for mutation in run_summary:
        if run_summary[mutation]['sample_count'] == num_samples and run_summary[mutation]['mixed_sample_count'] == 0:
            continue
        site_summary[mutation] = {
            'mutation':mutation,
            'mutation_freq':run_summary[mutation]['mutation_freq'],
            'sample_count':run_summary[mutation]['sample_count'],
            'mixed_sample_count':run_summary[mutation]['mixed_sample_count']
        }


    df = pd.DataFrame().from_dict(site_summary,orient='index')
    fig = px.bar(df,x='mutation',
                 y='mutation_freq',
                 color='sample_count',
                 title="Run Mutation Site Average Frequency",
                 labels={'mutation': 'Scheme mutation key',
                         'mutation_freq':'Average mutation kmer frequency',
                         'sample_count':'Number of samples with mutation site present'},
                 color_continuous_scale=px.colors.sequential.Viridis)
    fig.write_html(os.path.join(outdir,"kmer-coverage.run.html"))

    df = df.sort_values(by=['sample_count'],ascending=False)

    fig = px.bar(df,x='mutation',
                 y='sample_count',
                 color='mixed_sample_count',
                 title="Run Mutation Site Sample counts",
                 labels={'mutation': 'Scheme mutation key',
                         'sample_count':'Count of samples with feature present',
                         'mixed_sample_count':'Count of samples with mixed kmer states'},
                 color_continuous_scale=px.colors.sequential.Viridis)
    fig.write_html(os.path.join(outdir,"kmer-mixed.run.html"))

def calc_sample_kmer_distMat(kmer_results_df):
    return nan_compatible_kmer_pairwise_distmatrix(kmer_results_df)

def create_sample_comparison_plots(scheme_info,kmer_results_df,labels,mds_outfile,dendro_outfile,feature_outfile):
    plots = {
        'mds':'',
        'sample_dendro':'',
        'feature_dendro':'',
        'coverage':'',
        'clusters_assignments':{}
    }
    labels = kmer_results_df.columns.tolist()

    #Create kmer dist matrix
    disMat = calc_sample_kmer_distMat(kmer_results_df)

    # Create MDS plot
    plots['mds'] = plot_mds(disMat, labels,mds_outfile)

    #Create Dendrogram heatmap
    d = dendrogram_visualization()

    plots['sample_dendro'] =  d.build_dendrogram_heatmap(labels, kmer_results_df.T, dendro_outfile)
    #filter kmers to just those present in at least one sample and the alt state
    uid_to_state = scheme_info['uid_to_state']
    uid_to_dna_name = scheme_info['uid_to_dna_name']
    filter_list = []

    for uid in uid_to_state:
        if uid_to_state[uid] == 'alt':
            filter_list.append(uid)

    alt_kmer_result_df = kmer_results_df.iloc[filter_list, :]
    alt_kmer_result_df["sum"] = alt_kmer_result_df.sum(axis=1)
    alt_kmer_result_df = alt_kmer_result_df[alt_kmer_result_df["sum"] > 0]
    alt_kmer_result_df = alt_kmer_result_df.drop(['sum'], axis=1)

    dna_names = []


    # Create Dendrogram heatmap of features
    if len(alt_kmer_result_df) > 0:
        for uid in alt_kmer_result_df.index.values.tolist():
            dna_names.append("{}:{}".format(uid, uid_to_dna_name[uid]))
        plots['feature_dendro'] = d.build_dendrogram_heatmap_features(labels, dna_names, alt_kmer_result_df.T, feature_outfile)
    else:
        for uid in kmer_results_df.index.values.tolist():
            dna_names.append("{}:{}".format(uid, uid_to_dna_name[uid]))
        plots['feature_dendro'] = d.build_dendrogram_heatmap_features(labels, dna_names, kmer_results_df.T,
                                                                      feature_outfile)
    return plots



def run():
    cmd_args = parse_args()

    analysis_parameters = vars(cmd_args)

    logger = init_console_logger(2)
    is_args_ok = validate_args(cmd_args, logger)
    if not is_args_ok:
        logger.error("One or more command line arguments has an issue, please check the log messages and try again")
        sys.exit()

    # input parameters
    prefix = cmd_args.prefix
    min_cov = cmd_args.min_cov
    min_cov_frac = cmd_args.min_cov_frac
    max_missing_sites = cmd_args.max_missing_sites
    max_mixed_sites = cmd_args.max_mixed_sites
    scheme_file = cmd_args.scheme
    sample_type = cmd_args.type
    trim_seqs = cmd_args.trim
    R1 = cmd_args.R1
    R2 = cmd_args.R2
    SE = cmd_args.se
    seqTech = cmd_args.seqTech.lower()
    perform_read_correction = cmd_args.correct
    perform_read_dedup = cmd_args.dedup
    min_read_len = cmd_args.min_read_len
    no_plots = cmd_args.no_plots
    force_plots = cmd_args.force_plots
    data_dir = cmd_args.data_dir
    outdir = cmd_args.outdir
    n_threads = int(cmd_args.n_threads)
    primers = cmd_args.primers
    no_template_control = cmd_args.no_template_control
    positive_control = cmd_args.positive_control
    genotype_dist_cutoff = cmd_args.genotype_dist_cutoff
    min_genome_size = cmd_args.min_genome_size
    max_genome_size = cmd_args.max_genome_size
    min_genome_cov_depth = cmd_args.min_genome_cov_depth
    max_features = cmd_args.max_features
    type_only = cmd_args.typer_only
    scheme_meta = cmd_args.scheme_meta
    sample_meta = cmd_args.sample_meta

    # Initialize scheme
    if scheme_file in TYPING_SCHEMES:
        logger.info("User selected built-in scheme {}".format(scheme_file))
        scheme_name = scheme_file
        scheme_file = TYPING_SCHEMES[scheme_file]
    else:
        scheme_name = os.path.basename(scheme_file)

    logger.info("Reading kmer scheme from {}".format(scheme_file))
    scheme = parseScheme(scheme_file)
    logger.info("Initializing scheme data structure from {}".format(scheme_file))
    scheme_info = constructSchemeLookups(scheme)

    num_kmers = len(scheme_info['uid_to_kseq'])
    logger.info("Scheme interogates {} mutations".format(len(scheme)))
    logger.info("Scheme contains {} kmers".format(len(scheme_info['uid_to_kseq'])))
    logger.info("Min kmer len: {} , Max kmer len: {}".format(scheme_info['min_kmer_len'], scheme_info['max_kmer_len']))
    perform_scheme_check = False

    genotype_metadata = {}
    if scheme_meta is not None:
        scheme_meta_df = read_tsv(scheme_meta)
        columns = list(scheme_meta_df.columns)
        for row in scheme_meta_df.iterrows():
            genotype_metadata[row.genotype] = {}
            for field in columns:
                if field == 'genotype':
                    continue
                genotype_metadata[row.genotype][field] = row[field]

    sample_metadata = {}
    if sample_meta is not None:
        sample_meta_df = read_tsv(sample_meta)
        columns = list(sample_meta_df.columns)
        for row in sample_meta_df.iterrows():
            sample_metadata[row.sample_id] = {}
            for field in columns:
                if field == 'sample_id':
                    continue
                sample_metadata[row.sample_id][field] = row[field]

    ambiguousGenotypes = {}

    if perform_scheme_check:
        logger.info("Validating scheme genotyping ruleset")
        ambiguousGenotypes = detectAmbigGenotypes(scheme_info)

        if len(ambiguousGenotypes['kmers']) > 0:
            logger.warn("Scheme contains genotypes which cannot be reliably identified based on kmers")
            for geno in ambiguousGenotypes['kmers']:
                conflicts = [geno] + ambiguousGenotypes['kmers'][geno]
                conflicts = '; '.join([str(x) for x in conflicts])
                logger.warn("Ambiguous kmer genotypes {}".format(conflicts))

        if len(ambiguousGenotypes['mutations']) > 0:
            logger.warn("Scheme contains genotypes which cannot be reliably identified based on mutations")
            for geno in ambiguousGenotypes['mutations']:
                conflicts = [geno] + ambiguousGenotypes['mutations'][geno]
                conflicts = '; '.join([str(x) for x in conflicts])
                logger.warn("Ambiguous mutation genotypes {}".format(conflicts))

    # Identify if scheme kmers overlap any primer kmers
    primer_kmers = {}
    is_overlap = False
    overlapping_primers = {}
    if primers is not None:
        if primers in PRIMER_SCHEMES:
            primers = PRIMER_SCHEMES[primers]
        logger.info("Loading primer sequences from {}".format(primers))
        primer_df = read_tsv(primers)

        max_primer_len = 0
        for row in primer_df.itertuples():
            primer_kmers[row.name] = row.seq
            klen = len(row.seq)
            if klen > max_primer_len:
                max_primer_len = klen
        del (primer_df)
        logger.info("Determining longest primer length to use for read trimming: {}".format(max_primer_len))
        logger.info("Detecting kmers overlapping primers")
        overlapping_primers = checkPrimerKmerOverlap(scheme_info, primer_kmers)

        if len(overlapping_primers) > 0:
            is_overlap = True
        logger.info("There are {} kmers overlapping primers".format(len(overlapping_primers)))


    trim_front_bp = 0
    trim_tail_bp = 0
    if is_overlap and trim_seqs:
        trim_front_bp = max_primer_len
        trim_tail_bp = max_primer_len

    # initialize analysis directory
    if not os.path.isdir(outdir):
        logger.info("Creating analysis results directory {}".format(outdir))
        os.mkdir(outdir, 0o755)
    else:
        logger.info("Results directory {} already exits, will overwrite any results files here".format(outdir))

    # Logging
    analysis_date = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    # Gather sequence files in directory if specified
    seqManifest = {}
    if data_dir:
        logger.info("Identifying sequences in directory {} ".format(data_dir))
        seqManifest = create_seq_manifest(data_dir)
    elif R1 and R2:
        if R1[-2:] == 'gz':
            is_compressed = True
        seqManifest[prefix] = {
            'seq_type': 'fastq',
            'data_type': 'sample',
            'is_valid': True,
            'is_compressed': is_compressed,
            'is_paired': True,
            'seq_files': [R1, R2],
            'messages': []
        }
    elif SE:
        if SE[-2:] == 'gz':
            is_compressed = True
        seqManifest[prefix] = {
            'seq_type': 'fastq',
            'data_type': 'sample',
            'is_valid': True,
            'is_compressed': is_compressed,
            'is_paired': False,
            'seq_files': [SE],
            'messages': []
        }

    if no_template_control is not None:
        if no_template_control in seqManifest:
            seqManifest[no_template_control]['data_type'] = 'negative'

    if positive_control is not None:
        if positive_control in seqManifest:
            seqManifest[positive_control]['data_type'] = 'positive'

    if seqTech not in ['nanopore', 'illumina']:
        logging.error(
            "Please specify nanopore or illumina as the sequencing technology, mixed runs are not supported")

    if len(seqManifest) == 0:
        logger.error("Error no fasta/fastq files to process")
        sys.exit()

    sampleManifest = init_sample_manifest(seqManifest, scheme_info, scheme_name, analysis_date, sample_type, seqTech)

    # process reads
    if not type_only:
        read_dir = os.path.join(outdir, "processed_reads")
        if not os.path.isdir(read_dir):
            logger.info("Creating processed read directory {}".format(read_dir))
            os.mkdir(read_dir, 0o755)
        else:
            logger.info("Results directory {} already exits, will overwrite any results files here".format(read_dir))

        logger.info("Determine seq stats for {} samples".format(len(sampleManifest)))
        sampleManifest = process_reads(sampleManifest, perform_read_correction, read_dir, seqTech, read_dir,
                                       min_read_len, trim_front_bp, trim_tail_bp, perform_read_dedup,
                                       n_threads=n_threads)

    logger.info("Performing kmer detection for {} samples".format(len(sampleManifest)))
    kmer_counts = perform_kmer_searching(sampleManifest, scheme_info, min_cov, n_threads=n_threads)
    logger.info("Filtering low freq kmers for {} samples".format(len(sampleManifest)))
    kmer_counts = filter_low_freq_kmers(kmer_counts,num_kmers,min_freq=min_cov)
    logger.info("Filtering low frac kmers for {} samples".format(len(sampleManifest)))
    kmer_counts = filter_low_frac_kmers(kmer_counts,scheme_info,min_cov_frac)

    if no_template_control is not None:
        kmer_counts = filter_contam_kmers(kmer_counts, no_template_control, scheme_info, min_cov_frac)
    logger.info("Calculating genotype distances for {} samples accross {} genotypes".format(len(sampleManifest),len(scheme_info['genotypes'])))


    kmer_geno_dists = comp_kmer_profiles(kmer_counts, scheme_info)
    mutation_fracs = calc_site_frac(kmer_counts, scheme_info)

    mutation_geno_dists = comp_mutation_profiles(mutation_fracs, scheme_info)

    for sample_id in mutation_geno_dists:
        mutation_geno_dists[sample_id] = {k: v for k, v in sorted(mutation_geno_dists[sample_id].items(), key=lambda item: item[1])}

    for sample_id in kmer_geno_dists:
        kmer_geno_dists[sample_id] = {k: v for k, v in sorted(kmer_geno_dists[sample_id].items(), key=lambda item: item[1])}


    logger.info("Assigning samples to geneotypes")
    genotype_assignments = assign_genotypes(kmer_counts, mutation_geno_dists, kmer_geno_dists, scheme_info, genotype_dist_cutoff)

    logging.info("Detecting mixed sites")
    mixed_sites = get_mixed_sites(mutation_fracs, scheme_info,min_cov_frac)
    for sample_id in mixed_sites:
        sampleManifest[sample_id]['num_detected_scheme_mixed_mutations'] = len(mixed_sites[sample_id])
        sampleManifest[sample_id]['detected_scheme_mixed_mutations'] = mixed_sites[sample_id]

    if not type_only:
        logging.info("Calculating genome size for samples")
        sampleManifest = calc_genome_size(sampleManifest, outdir, min_genome_cov_depth, kLen=21, n_threads=n_threads)

    logging.info("Adding kmer results to manifest")
    sampleManifest = add_kmer_count_data(sampleManifest, mutation_fracs, scheme_info, kmer_counts)

    logging.info("Adding genotyping results to manifest")
    sampleManifest = add_genotype_info(sampleManifest, scheme_info, mutation_fracs, genotype_assignments)

    logging.info("Performing QA on samples")
    sampleManifest = QA_results(sampleManifest, min_genome_cov_depth, max_missing_sites, min_genome_size, max_genome_size,max_mixed_sites)

    # create sample summary table
    write_sample_summary_results(sampleManifest, scheme_info, os.path.join(outdir,"{}.sample.composition.report.txt".format(prefix)), sample_metadata,genotype_metadata, max_features)

    #write kmer profile
    kmer_df = pd.DataFrame(kmer_counts)
    kmer_df.reset_index().to_csv(os.path.join(outdir,"{}-kmer.profile".format(prefix)),index=False,sep="\t")

    if not type_only and len(sampleManifest) > 1:
        # create run summary table
        run_summary = create_run_summary(sampleManifest, kmer_counts, scheme_info)

        # write run summary
        write_run_summary(run_summary, outdir, len(sampleManifest), prefix)

        labels = kmer_df.columns.tolist()
        if len(labels) > 1:
            mds_outfile = os.path.join(outdir, "{}-mds.html".format(prefix))
            dendro_outfile = os.path.join(outdir, "{}-sample.dendrogram.html".format(prefix))
            feature_outfile = os.path.join(outdir, "{}-feature.heatmap.html".format(prefix))
            create_sample_comparison_plots(scheme_info, kmer_df, labels, mds_outfile, dendro_outfile, feature_outfile)

