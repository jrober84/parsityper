#!/usr/bin/python
import glob
import time
from argparse import (ArgumentParser, FileType)
import logging, os, sys, re, collections, operator, math, shutil, datetime
from collections import Counter
import pandas as pd
from parsityper.helpers import validate_args, init_console_logger, read_tsv, scheme_to_biohansel_fasta, process_biohansel_kmer, \
    init_kmer_targets, get_scheme_template, generate_random_phrase, calc_md5, \
    get_kmer_groups, get_kmer_group_mapping, get_sequence_files, read_genotype_profiles, dist_compatible_profiles, \
    generate_target_presence_table, nan_compatible_kmer_pairwise_distmatrix
import copy
from parsityper.bio_hansel import bio_hansel
import statistics
from parsityper.constants import PRIMER_SCHEMES, TYPING_SCHEMES
from parsityper.helpers import profile_pairwise_distmatrix
from parsityper.visualizations import dendrogram_visualization, plot_mds, generate_sample_coverage_plot
from parsityper.db_search import process_strain_db, get_neighbours

from parsityper.ext_tools.kmc import kmc_summary
from parsityper.ext_tools.canu import run_canu_correct
from parsityper.ext_tools.lighter import run_lighter
from parsityper.ext_tools.fastp import run_fastp
from parsityper.ext_tools.mash import mash_sample_comparisons
from parsityper.scheme import parseScheme, constructSchemeLookups, detectAmbigGenotypes
from parsityper.kmerSearch.kmerSearch import init_automaton_dict,perform_kmerSearch_fasta,perform_kmerSearch_fastq,process_kmer_results
from parsityper.helpers import read_fasta
from multiprocessing import Pool

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Parsimony based sample Kmer-analysis')
    parser.add_argument('--data_dir', type=str, required=False,
                        help='directory of fasta/fastq files')
    parser.add_argument('--outdir', type=str, required=True,
                        help='output directory')
    parser.add_argument('--prefix', type=str, required=False,
                        help='output file prefix', default='parsityper_analysis')
    parser.add_argument('--scheme', type=str, required=True,
                        help='TSV formated kmer scheme', default='')
    parser.add_argument('--primers', type=str, required=False,
                        help='TSV formated primer file', default=None)
    parser.add_argument('--mode', type=str, required=True,
                        help='Operate in batch or single mode. Incompatible with input options other than data_dir')
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
    parser.add_argument('--se', type=str, required=False,
                        help='single-end fastq read')
    parser.add_argument('--R1', type=str, required=False,
                        help='paired-end fwd fastq read')
    parser.add_argument('--R2', type=str, required=False,
                        help='paired-end rev fastq read')
    parser.add_argument('--min_read_len', type=int, required=False,
                        help='Filter reads shorter than this limit',
                        default=0)
    parser.add_argument('--max_read_len', type=int, required=False,
                        help='Filter reads longer than this limit',
                        default=0)
    parser.add_argument('--min_cov', type=int, required=False,
                        help='Absolute minimum kmer coverage for fastq read detection default=auto determine coverage',
                        default=50)
    parser.add_argument('--min_cov_frac', type=float, required=False,
                        help='Minimum percentage of total pool required for k-mer detection range 0 - 1.0 (default=0.05)',
                        default=0.05)
    parser.add_argument('--max_mixed_sites', type=int, required=False,
                        help='Maximum number of sites allowed to have both kmer states', default=10)
    parser.add_argument('--max_missing_sites', type=int, required=False,
                        help='Maximum number of sites allowed to be missing', default=500)
    parser.add_argument('--sample_dist_cutoff', type=float, required=False,
                        help='Maximum profile distance for including sample in a comparison', default=0.05)
    parser.add_argument('--genotype_dist_cutoff', type=float, required=False,
                        help='Maximum profile distance for considering a genotype result valid', default=0.05)
    parser.add_argument('--n_threads', type=str, required=False,
                        help='output directory', default=1)
    parser.add_argument('--strain_profiles', type=str, required=False,
                        help='Kmer profile scheme to compare against')
    parser.add_argument('--no_template_control', type=str, required=False,
                        help='SampleID of the negative control sample')
    parser.add_argument('--positive_control', type=str, required=False,
                        help='Fasta/Fastq formatted positive control data')
    parser.add_argument('--list', type=str, required=False,
                        help='list in-built primer and typing schemes')
    parser.add_argument('--no_plots', required=False,
                        help='suppress making plots, required for large datasets',action='store_true')
    parser.add_argument('--force_plots', required=False,
                        help='Try plotting for datasets > 1000',action='store_true')
    parser.add_argument('--delete_processed_reads', required=False,
                        help='Delete processed reads after completion',action='store_true')
    return parser.parse_args()

def find_seq_files(input_dir):
    file_dict = {}
    fasta_file_extensions = ['.fasta','.fas','.fa','.fna']
    fastq_file_extensions = ['.fastq','.fq']
    paired_end_notations = ['_R1','_R2','_1','_2',]
    file_list = glob.glob("{}*".format(input_dir))
    for file in file_list:
        fileType = ''
        fileName = os.path.basename(file)
        sampleName = fileName
        for ext in fasta_file_extensions:
            if ext in file:
                fileType = 'fasta'
                sampleName = re.sub(r"{}$".format(ext), '', sampleName )
                break

        if fileType == '':
            for ext in fastq_file_extensions:
                if ext in file:
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

def compare_sample_to_genotypes(data,genotypes,kmer_profiles,genotype_pos_kmers,genotype_inf_kmers,genotype_par_kmers,mutation_to_uids,uid_to_mutation,uid_to_state,num_kmers,min_freq,min_cov_frac):
    genotype_results = {}
    for genotype in genotypes:
        genotype_results[genotype] = {
            'scheme_pos_kmers': len(genotype_pos_kmers[genotype]),
            'scheme_inf_kmers': len(genotype_inf_kmers[genotype]),
            'scheme_par_kmers': len(genotype_par_kmers[genotype]),
            'kmer_genotype_dist': 0,
            'num_informative_match': 0,
            'num_informative_missing': 0,
            'num_pos_match': 0,
            'num_pos_mismatch': 0,
            'matched_pos_kmers': [],
            'mismatched_pos_kmers': [],
            'mismatched_kmers': [],
            'is_compatible': True
        }

    for genotype in genotypes:
        dist = 0
        for mutation_key in mutation_to_uids:
            k_ids = mutation_to_uids[mutation_key]
            count_ref = 0
            count_alt = 0
            for k in k_ids:
                kFreq = data['raw_kmer_freq'][k]
                state = uid_to_state[k]
                if kFreq < min_freq:
                    continue
                if state == 'ref':
                    count_ref += kFreq
                else:
                    count_alt += kFreq
            total = count_ref + count_alt
            if total < min_freq:
                continue

            for k in k_ids:
                kFreq = data['raw_kmer_freq'][k]
                if kFreq < min_freq:
                    continue
                frac = kFreq / total

                if frac < min_cov_frac:
                    char = 0
                elif frac >= min_cov_frac and frac <= 1 - min_cov_frac:
                    char = 0.5
                else:
                    char = 1
                pChar = kmer_profiles[genotype][k]

                if char == pChar:
                    if char == 1:
                        genotype_results[genotype]['num_pos_match'] += 1
                        genotype_results[genotype]['matched_pos_kmers'].append(k)
                else:
                    if (char != 0.5 and pChar != 0.5) and (char != pChar):
                        dist += 1
                        genotype_results[genotype]['mismatched_kmers'].append(k)
                        if pChar == 1 and char == 0:
                            genotype_results[genotype]['num_pos_mismatch'] += 1
                            genotype_results[genotype]['mismatched_pos_kmers'].append(k)
                            genotype_results[genotype]['is_compatible'] = False
                        if pChar == 0 and char == 1:
                            genotype_results[genotype]['is_compatible'] = False

        genotype_results[genotype]['kmer_genotype_dist'] = dist / num_kmers
        genotype_results[genotype]['num_pos_mismatch'] = len(
        genotype_results[genotype]['mismatched_pos_kmers'])
        print("{}\t{}".format(genotype,genotype_results[genotype]['is_compatible']))

    return genotype_results

def summarize_genotype_kmers(scheme_info,kmer_results,min_freq,min_cov_frac,n_threads=1):
    kmer_profiles = scheme_info['kmer_profiles']
    genotype_pos_kmers = {}
    genotype_inf_kmers = {}
    genotype_par_kmers = {}

    for genotype in kmer_profiles:
        genotype_pos_kmers[genotype] = []
        genotype_inf_kmers[genotype] = []
        genotype_par_kmers[genotype] = []
        profile = kmer_profiles[genotype]
        for i in range(0,len(profile)):
            value = profile[i]
            if value == 1:
                genotype_pos_kmers[genotype].append(i)
            if value != 0.5:
                genotype_inf_kmers[genotype].append(i)
            else:
                genotype_inf_kmers[genotype].append(i)

    #initialize data
    genotypes = list(kmer_profiles.keys())
    samples = list(kmer_results.keys())
    genotype_results = {}
    for sample_id in samples:
        genotype_results[sample_id] = {}
        for genotype in genotypes:
            genotype_results[sample_id][genotype] = {
                'scheme_pos_kmers':len(genotype_pos_kmers[genotype]),
                'scheme_inf_kmers': len(genotype_inf_kmers[genotype]),
                'scheme_par_kmers': len(genotype_par_kmers[genotype]),
                'kmer_genotype_dist':0,
                'num_informative_match': 0,
                'num_informative_missing': 0,
                'num_pos_match': 0,
                'num_pos_mismatch': 0,
                'matched_pos_kmers': [],
                'mismatched_pos_kmers': [],
                'mismatched_kmers':[],
                'is_compatible':True
            }

    mutation_to_uids = scheme_info['mutation_to_uid']
    uid_to_mutation = scheme_info['uid_to_mutation']
    uid_to_state =  scheme_info['uid_to_state']
    num_kmers = len(uid_to_mutation)

    if n_threads > 1:
        pool = Pool(processes=n_threads)

    for sample_id in samples:
        data = kmer_results[sample_id]
        if n_threads == 1:
            genotype_results[sample_id] = compare_sample_to_genotypes(data,genotypes,kmer_profiles,
                                                                      genotype_pos_kmers,genotype_inf_kmers,
                                                                      genotype_par_kmers,mutation_to_uids,
                                                                      uid_to_mutation,uid_to_state,num_kmers,
                                                                      min_freq,min_cov_frac)
        else:
            genotype_results[sample_id] = pool.apply_async(compare_sample_to_genotypes,
                                                           (data,genotypes,kmer_profiles,
                                                                      genotype_pos_kmers,genotype_inf_kmers,
                                                                      genotype_par_kmers,mutation_to_uids,
                                                                      uid_to_mutation,uid_to_state,num_kmers,
                                                                      min_freq,min_cov_frac))

    if n_threads > 1:
        pool.close()
        pool.join()
        for sample_id in genotype_results:
            genotype_results[sample_id] = genotype_results[sample_id].get()
    return genotype_results


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

def call_compatible_genotypes(scheme_info, kmer_results, min_freq, min_cov_frac,n_threads=1):
    logging.info("Comparing kmer profiles against known genotypes")
    genotype_results = summarize_genotype_kmers(scheme_info, kmer_results, min_freq, min_cov_frac,n_threads)
    called_genotypes = {}
    logging.info("Identifying compatible genotypes")
    for sample_id in genotype_results:
        valid_genotypes = {}
        parsiony_genotypes = {}
        called_genotypes[sample_id] = {'called_genotypes':{},'all_genotypes':genotype_results[sample_id]}
        for genotype in genotype_results[sample_id]:
            if not genotype_results[sample_id][genotype]:
                continue
            valid_genotypes[genotype] = genotype_results[sample_id][genotype]['kmer_genotype_dist']

        #order valid genotypes by lowest kmer distance
        valid_genotypes = {k: v for k, v in sorted(valid_genotypes.items(), key=lambda item: item[1])}
        assigned_kmers = []
        for genotype in valid_genotypes:
            unassigned = list(set(genotype_results[sample_id][genotype]['matched_pos_kmers'] )- set(assigned_kmers))
            if len(unassigned) == 0:
                continue
            assigned_kmers.extend(unassigned)
            parsiony_genotypes[genotype] = unassigned
        called_genotypes[sample_id]['called_genotypes'] = parsiony_genotypes

    return called_genotypes

def run():
    cmd_args = parse_args()
    logger = init_console_logger(2)
    is_args_ok = validate_args(cmd_args, logger)
    if not is_args_ok:
        logger.error("One or more command line arguments has an issue, please check the log messages and try again")
        sys.exit()

    # input parameters
    prefix = cmd_args.prefix
    min_cov = cmd_args.min_cov
    min_cov_frac = cmd_args.min_cov_frac
    detection_limit = min_cov
    max_missing_sites = cmd_args.max_missing_sites
    max_mixed_sites = cmd_args.max_mixed_sites
    scheme_file = cmd_args.scheme
    mode = cmd_args.mode
    type = cmd_args.type
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
    nthreads = int(cmd_args.n_threads)
    strain_profiles = cmd_args.strain_profiles
    primers = cmd_args.primers
    no_template_control = cmd_args.no_template_control
    positive_control = cmd_args.positive_control
    genotype_dist_cutoff = cmd_args.genotype_dist_cutoff


    #Initialize scheme
    logger.info("Reading kmer scheme from {}".format(scheme_file))
    scheme = parseScheme(scheme_file)
    logger.info("Initializing scheme data structure from {}".format(scheme_file))
    scheme_info = constructSchemeLookups(scheme)
    logger.info("Scheme interogates {} mutations".format(len(scheme)))
    logger.info("Scheme contains {} kmers".format(len(scheme_info['uid_to_kseq'])))
    logger.info("Min kmer len: {} , Max kmer len: {}".format(scheme_info['min_kmer_len'],scheme_info['max_kmer_len']))
    perform_scheme_check = False
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

    #Init Ahocorasak automation objects
    logger.info("Initializing aho-corasick automation")
    aho = {'scheme':init_automaton_dict(scheme_info['uid_to_kseq'])}

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
        del(primer_df)
        logger.info("Determining longest primer length to use for read trimming: {}".format(max_primer_len))
        logger.info("Detecting kmers overlapping primers")
        overlapping_primers = checkPrimerKmerOverlap(scheme_info, primer_kmers)

        if len(overlapping_primers) > 0:
            is_overlap = True
        logger.info("There are {} kmers overlapping primers".format(len(overlapping_primers)))
        aho['primers'] = init_automaton_dict(primer_kmers)

    if is_overlap:
        trim_front_bp = max_primer_len
        trim_end_bp = max_primer_len

    # initialize analysis directory
    if not os.path.isdir(outdir):
        logger.info("Creating analysis results directory {}".format(outdir))
        os.mkdir(outdir, 0o755)
    else:
        logger.info("Results directory {} already exits, will overwrite any results files here".format(outdir))

    # output filenames
    report_run_info_log = open(os.path.join(outdir, "{}.run.info.log".format(prefix)), 'w')
    report_sample_composition_detailed = os.path.join(outdir,
                                                      "{}.sample_composition.report.detailed.txt".format(prefix))
    report_sample_composition_summary = os.path.join(outdir, "{}.sample_composition.report.summary.txt".format(prefix))
    report_positive = os.path.join(outdir, "{}.positive_control.report.txt".format(prefix))
    report_negative = os.path.join(outdir, "{}.negative_control.report.txt".format(prefix))
    report_primers = os.path.join(outdir, "{}.primers.report.txt".format(prefix))
    report_run_info_composition = os.path.join(outdir, "{}.run.info.report.txt".format(prefix))
    report_run_kmer_composition = os.path.join(outdir, "{}.run_composition.report.txt".format(prefix))
    report_qc_metrics = os.path.join(outdir, "{}.run.qc.txt".format(prefix))
    report_individual_sample_html_prefix = os.path.join(outdir, "{}.sample_#.report.html".format(prefix))
    report_genotype_targets = os.path.join(outdir, "{}.sample.genotype.targets.txt".format(prefix))
    report_dbsearch = os.path.join(outdir, "{}.dbsearch.txt".format(prefix))

    #Images
    report_run_kmer_dendrogram = os.path.join(outdir, "{}.run_kmer_dendrogram.png".format(prefix))
    report_scipy_profile = os.path.join(outdir, "{}.scipy.features.profile".format(prefix))
    report_run_kmer_mds = os.path.join(outdir, "{}.run_kmer_mds.png".format(prefix))
    report_run_kmer_coverage = os.path.join(outdir, "{}.run_kmer_coverage.png".format(prefix))

    report_run_info_log.write("Start Time\t{}\n".format(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
    report_run_info_log.write("Output Directory\t{}\n".format(outdir))
    report_run_info_log.write("Mode\t{}\n".format(mode))
    report_run_info_log.write("Type\t{}\n".format(type))
    report_run_info_log.write("Distance Cutoff\t{}\n".format(genotype_dist_cutoff))
    report_run_info_log.write("Min kmer freq\t{}\n".format(min_cov))
    report_run_info_log.write("Min kmer ratio\t{}\n".format(min_cov_frac))
    report_run_info_log.write("Detection limit\t{}\n".format(detection_limit))
    report_run_info_log.write("Max missing sites\t{}\n".format(max_missing_sites))
    report_run_info_log.write("Max mixed sites\t{}\n".format(max_missing_sites))
    report_run_info_log.write("Strain Profiles\t{}\n".format(strain_profiles))
    report_run_info_log.write("Positive Control\t{}\n".format(positive_control))
    report_run_info_log.write("No Template Control\t{}\n".format(no_template_control))
    report_run_info_log.write("Primers\t{}\n".format(primers))

    # Perform read correction if requested
    # initialize read directory
    if perform_read_correction:
        read_dir = os.path.join(outdir, "processed_reads")
        if not os.path.isdir(read_dir):
            logger.info("Creating processed read directory {}".format(read_dir))
            os.mkdir(read_dir, 0o755)
        else:
            logger.info("Results directory {} already exits, will overwrite any results files here".format(read_dir))

    #Gather sequence files in directory if specified
    seqManifest = {}
    if data_dir:
        logger.info("Identifying sequences in directory {} ".format(data_dir))
        seqManifest = create_seq_manifest(data_dir)
    elif R1 and R2:
        if  R1[-2:] == 'gz':
            is_compressed = True
        seqManifest[prefix] =  {
            'seq_type':'fastq',
            'data_type': 'sample',
            'is_valid': True,
            'is_compressed': is_compressed,
            'is_paired': True,
            'seq_files': [R1,R2],
            'messages': []
        }
    elif SE:
        if  SE[-2:] == 'gz':
            is_compressed = True
        seqManifest[prefix] =  {
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
        logging.error("Please specify nanopore or illumina as the sequencing technology, mixed runs are not supported")

    #Process the identified reads and determine summary stats for the sample
    sampleManifest = {}

    for sampleID in seqManifest:
        sampleManifest[sampleID] = {
            'seq_type': seqManifest[sampleID]['seq_type'],
            'est_genome_size':0,
            'num_unique_kmers':0,
            'num_counted_unique_kmers':0,
            'raw_seq_files': seqManifest[sampleID]['seq_files'],
            'num_reads': len(seqManifest[sampleID]['seq_files']),
            'total_reads_pre':0,
            'total_reads_post': 0,
            'total_bases_pre': 0,
            'total_bases_post': 0,
            'read_mean_len_pre':0,
            'read_mean_len_post': 0,
            'gc_pre': 0,
            'gc_post': 0,
            'insert_size_peak': 0,
            'duplication_rate': 0,
            'processed_reads':[]
        }
        if sampleManifest[sampleID]['seq_type'] == 'fasta':
            continue
        # get GenomeSize using KMC based on fwd read
        if sampleManifest[sampleID]['seq_type'] == 'fastq':
            data = kmc_summary(seqManifest[sampleID]['seq_files'][0], os.path.join(outdir,"kmc.tmp1"), os.path.join(outdir,"kmc.tmp2"), freq=10, kmer_length=21, n_threads=1)
            os.remove(os.path.join(outdir,"kmc.tmp1"))
            os.remove(os.path.join(outdir, "kmc.tmp2"))
            genomeSize = data['est_genome_size']

        for field in data:
            sampleManifest[sampleID][field] = data[field]

        # Perform read correction if requested
        if perform_read_correction:
            # initialize read directory
            sample_read_dir = os.path.join(read_dir,"{}".format(sampleID))
            if not os.path.isdir(sample_read_dir):
                logger.info("Creating processed read directory {}".format(sample_read_dir))
                os.mkdir(sample_read_dir, 0o755)
            else:
                logger.info("Results directory {} already exits, will overwrite any results files here".format(sample_read_dir))

            if seqTech == 'illumina':
                run_lighter(seqManifest[sampleID]['seq_files'],sample_read_dir,17,genomeSize,n_threads=nthreads)
                sampleManifest[sampleID]['processed_reads'] = glob.glob("{}/*.cor.*".format(sample_read_dir))

            else:
                run_canu_correct(seqManifest[sampleID]['seq_files'][0], sampleID, sample_read_dir, genomeSize, min_length=500, minOverlapLength=500,
                                 corOutCoverage=1000,
                                 n_threads=nthreads)
                #TODO add the canu specific corrected reads
                sampleManifest[sampleID]['processed_reads'] = []

        #Perform read preprocessing using fastp
        if len(sampleManifest[sampleID]['processed_reads'])> 0:
            read_set = sampleManifest[sampleID]['processed_reads']
        else:
            read_set = sampleManifest[sampleID]['raw_seq_files']

        fastp_dir = os.path.join(outdir,"fastp")
        if not os.path.isdir(sample_read_dir):
            logger.info("Creating processed read directory {}".format(sample_read_dir))
            os.mkdir(sample_read_dir, 0o755)
        else:
            logger.info(
                "Results directory {} already exits, will overwrite any results files here".format(sample_read_dir))
        merge = False
        if len(read_set) == 2:
            merge = True

        fastp_results = run_fastp(read_set,  fastp_dir, sampleID, min_read_len=min_read_len, trim_front_bp=trim_front_bp, trim_tail_bp=trim_end_bp, report_only=False,
                  dedup=False, merge_reads=merge, n_threads=nthreads)

        sampleManifest[sampleID]['total_reads_pre'] = fastp_results['summary']["before_filtering"]["total_reads"]
        sampleManifest[sampleID]['total_bases_pre'] = fastp_results['summary']["before_filtering"]["total_bases"]
        sampleManifest[sampleID]['gc_pre'] = fastp_results['summary']["before_filtering"]["gc_content"]
        sampleManifest[sampleID]['total_reads_post'] = fastp_results['summary']["after_filtering"]["total_reads"]
        sampleManifest[sampleID]['total_bases_post'] = fastp_results['summary']["after_filtering"]["total_bases"]
        sampleManifest[sampleID]['gc_post'] = fastp_results['summary']["after_filtering"]["gc_content"]
        sampleManifest[sampleID]['read_mean_len_pre'] = fastp_results['summary']["before_filtering"]["read1_mean_length"]
        sampleManifest[sampleID]['read_mean_len_post'] = fastp_results['summary']["after_filtering"]["read1_mean_length"]
        sampleManifest[sampleID]['duplication_rate'] = fastp_results["duplication"]

        if len(read_set) == 2:
            sampleManifest[sampleID]['read_mean_len_pre'] = (fastp_results[sampleID]['read_mean_len_pre'] +
                                                            fastp_results['summary']["before_filtering"][
                                                                "read2_mean_length"]) / 2
            sampleManifest[sampleID]['read_mean_len_post'] = (fastp_results[sampleID]['read_mean_len_post'] +
                                                            fastp_results['summary']["after_filtering"][
                                                                "read2_mean_length"]) / 2

            pR1 = os.path.join(fastp_dir,"{}_1.fastq".format(sampleID))
            pR2 = os.path.join(fastp_dir, "{}_2.fastq".format(sampleID))
            pM = os.path.join(fastp_dir,"{}.merge.fastq".format(sampleID))
            sampleManifest[sampleID]['processed_reads'] = [pR1,pR2,pM]

    kmer_results = {}
    #Identify kmers in each sample
    if nthreads > 1:
        pool = Pool(processes=nthreads)
    logger.info("Performing kmer searching on {} samples with {} threads".format(len(sampleManifest),nthreads))

    for sampleID in sampleManifest:
        if len(sampleManifest[sampleID]['processed_reads'])> 0:
            read_set = sampleManifest[sampleID]['processed_reads']
        else:
            read_set = sampleManifest[sampleID]['raw_seq_files']

        if nthreads == 1:
            if sampleManifest[sampleID]['seq_type'] == 'fastq':
                kmer_results[sampleID] = perform_kmerSearch_fastq(scheme_info['uid_to_kseq'], scheme_info['kseq_to_uids'],aho['scheme'], read_set)
            else:
                kmer_results[sampleID] = perform_kmerSearch_fasta(scheme_info['uid_to_kseq'], scheme_info['kseq_to_uids'],aho['scheme'], read_fasta(read_set[0]), min_cov)
        else:
            if sampleManifest[sampleID]['seq_type'] == 'fastq':
                kmer_results[sampleID] = pool.apply_async(perform_kmerSearch_fastq, (scheme_info['uid_to_kseq'], scheme_info['kseq_to_uids'],aho['scheme'], read_set))
            else:
                kmer_results[sampleID] = pool.apply_async(perform_kmerSearch_fasta,
                                                          (scheme_info['uid_to_kseq'], scheme_info['kseq_to_uids'],aho['scheme'], read_fasta(read_set[0]), min_cov))


    #Extract results in multithreaded mode
    if nthreads > 1:
        pool.close()
        pool.join()
        for sampleID in sampleManifest:
            kmer_results[sampleID] = kmer_results[sampleID].get()

    logger.info("Kmer searching complete")
    logger.info("Processing kmer results")

    #process kmer results
    kmer_results = process_kmer_results(scheme_info, kmer_results, min_cov, min_cov_frac)

    #Compare kmer results with genotype rules in the scheme
    logger.info("Comparing kmer results to genotypes")
    sample_genotype_results = call_compatible_genotypes(scheme_info, kmer_results, min_cov, min_cov_frac,nthreads)

    for sample_id in sample_genotype_results:
        print("{}\t{}".format(sample_id,sample_genotype_results[sample_id]['called_genotypes'].keys()))


run()




