#!/usr/bin/python
import glob
from parsityper.version import __version__
from argparse import (ArgumentParser, FileType)
import logging, os, sys, re, datetime, statistics
from pathlib import Path
import pandas as pd
from parsityper.helpers import validate_args, init_console_logger, read_tsv, \
generate_random_phrase, calc_md5, nan_compatible_kmer_pairwise_distmatrix
import copy
from parsityper.constants import PRIMER_SCHEMES, TYPING_SCHEMES, TYPER_SAMPLE_SUMMARY_HEADER_BASE, FIGURE_CAPTIONS, BATCH_HTML_REPORT
from parsityper.visualizations import dendrogram_visualization, plot_mds, create_heatmap
from jinja2 import Template
from parsityper.ext_tools.kmc import kmc_summary
from parsityper.ext_tools.canu import run_canu_correct
from parsityper.ext_tools.lighter import run_lighter
from parsityper.ext_tools.fastp import run_fastp
from parsityper.scheme import parseScheme, constructSchemeLookups, detectAmbigGenotypes
from parsityper.kmerSearch.kmerSearch import init_automaton_dict,perform_kmerSearch_fasta,perform_kmerSearch_fastq,process_kmer_results
from parsityper.helpers import read_fasta
from multiprocessing import Pool
import plotly.express as px
from plotly.offline import plot

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Parsimony based sample Kmer-analysis')
    parser.add_argument('--data_dir', type=str, required=False,
                        help='directory of fasta/fastq files')
    parser.add_argument('--sample_manifest', type=str, required=False,
                        help='TSV file of samples and files (sample_id,seq_file1,seq_file_2')
    parser.add_argument('--outdir', type=str, required=True,
                        help='output directory')
    parser.add_argument('--prefix', type=str, required=False,
                        help='output file prefix', default='parsityper_analysis')
    parser.add_argument('--scheme', type=str, required=True,
                        help='TSV formated kmer scheme', default='')
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
    parser.add_argument('--min_genome_size', type=int, required=False,
                        help='Flag genomes which are too small',
                        default=0)
    parser.add_argument('--max_genome_size', type=int, required=False,
                        help='Flag genomes which are too large',
                        default=-1)
    parser.add_argument('--min_genome_cov_depth', type=int, required=False,
                        help='Flag samples where the raw sequencing depth is too low',
                        default=40)
    parser.add_argument('--min_cov', type=int, required=False,
                        help='Absolute minimum kmer coverage for fastq read detection default=auto determine coverage',
                        default=50)
    parser.add_argument('--min_cov_frac', type=float, required=False,
                        help='Minimum percentage of total pool required for k-mer detection range 0 - 1.0 (default=0.05)',
                        default=0.05)
    parser.add_argument('--max_mixed_sites', type=int, required=False,
                        help='Maximum number of sites allowed to have both kmer states', default=100)
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
                        help='SampleID of the positive control sample')
    parser.add_argument('--list', type=str, required=False,
                        help='list in-built primer and typing schemes')
    parser.add_argument('--no_plots', required=False,
                        help='suppress making plots, required for large datasets',action='store_true')
    parser.add_argument('--force_plots', required=False,
                        help='Try plotting for datasets > 1000',action='store_true')
    parser.add_argument('--delete_processed_reads', required=False,
                        help='Delete processed reads after completion',action='store_true')
    parser.add_argument('--typer_only', required=False,
                        help='Skip read preprocessing and stat calcuations', action='store_true')
    parser.add_argument('--max_features', type=str, required=False,
                        help='max gene features to report', default=15)
    parser.add_argument('-V', '--version', action='version', version='%(prog)s {}'.format(__version__))

    return parser.parse_args()

def batch_sample_html_report(htmlStr,data,outfile):
    t = Template(htmlStr)
    html = t.render(data)
    # to save the results
    with open(outfile, "w") as fh:
        fh.write(html )
    fh.close()

def prep_hmtml_report_data(plots,sample_data):
    report = {
        'analysis_date':'',
        'sample_type':'',
        'num_strain_detected':'',
        'ave_strain_kmer_freq':'',
        'ave_mixed_kmer_freq':'',
        'total_samples':'',
        'count_failed_samples':'',
        'failed_sample_ids':'',
        'coverage_plot':'',
        'coverage_plot_caption':'',
        'mixed_target_plot':'',
        'mixed_target_plot_caption':'',
        'genotypes_found':'',
        'num_kmer_profiles':'',
        'genotype_abundance_plot':'',
        'genotype_abundance_plot_caption': '',
        'num_mutation_detected':'',
        'mutations_found':'',
        'mutation_table':'',
        'analysis_parameters':{}

    }
    return

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

def compare_sample_to_genotypes(data,genotype_results,scheme_info,outfile):
    uids = set(list(scheme_info['uid_to_mutation'].keys()))
    valid_uids = set(data['valid_uids'])
    mixed_sites = list(data['mixed_sites'])
    missing_sites = list(uids - valid_uids )
    exclude_sites = set(missing_sites + mixed_sites)
    mixed_sites = set(mixed_sites)
    missing_sites = set(missing_sites)

    detected_scheme_kmers = set(data['detected_scheme_kmers'])
    geno_rules = scheme_info['genotype_rule_sets']
    genotypes = list(genotype_results.keys())

    valid_genotypes = {}

    for mutation_key in scheme_info['mutation_to_uid']:
        if mutation_key in mixed_sites:
            continue
        detected_mut_kmers = scheme_info['mutation_to_uid'][mutation_key] & detected_scheme_kmers
        if len(detected_mut_kmers) == 0:
            continue
        uids = scheme_info['mutation_to_uid'][mutation_key]
        missing_mut_kmer = list(uids - detected_mut_kmers)
        detected_mut_kmers = list(detected_mut_kmers)
        state = scheme_info['uid_to_state'][detected_mut_kmers[0]]
        uids_to_add = []
        for uid in missing_mut_kmer:
            if state != scheme_info['uid_to_state'][uid]:
                continue
            uids_to_add.append(uid)
        exclude_sites = exclude_sites.union(set(uids_to_add))
    detected_scheme_kmers = detected_scheme_kmers - exclude_sites
    for genotype in genotypes:
        dist = 0
        informative_uids = set(geno_rules[genotype]['positive_uids']) | set(geno_rules[genotype]['partial_uids']) - exclude_sites
        matched = list(detected_scheme_kmers & informative_uids )
       # atypical = list(detected_scheme_kmers - informative_uids)
       # filt_atypical = []
       # for uid in atypical:
       #     mutation_key = scheme_info['uid_to_mutation'][uid]
       #     uids = set(scheme_info['mutation_to_uid'][mutation_key])
        #    ovl = set(geno_rules[genotype]['positive_uids'])
        #    if len(ovl) == 0:
        #        continue
        #    g_state = scheme_info['mutation_positions'][mutation_key]
        #    a_state = scheme_info['uid_to_state'][uid]
        #    if g_state == 0.5:
        #        continue
        #    if g_state == 0.5 or (g_state == 1 and a_state == 'alt') or (g_state == 0 and a_state == 'ref'):
        #        continue
        #    filt_atypical.append(uid)

       # atypical = filt_atypical
        mismatches = list(set(geno_rules[genotype]['positive_uids']) - set(matched) - exclude_sites)
        filtered_mismatches = []
        for uid in mismatches:
            state = scheme_info['uid_to_state'][uid]
            if uid in geno_rules[genotype]['positive_uids'] :
                filtered_mismatches.append(uid)


        mismatches = filtered_mismatches #+ filt_atypical
        #uids = (valid_uids & set(geno_rules[genotype]['positive_uids'])) - exclude_sites
        #uids = uids - set(geno_rules[genotype]['partial_uids'])
        #mismatches = uids - detected_scheme_kmers
        #matched = list(detected_scheme_kmers & (uids - mismatches))
        positive_alt_match = list(set(geno_rules[genotype]['positive_alt']) & set(matched ))
        positive_alt_mismatch = list(set(geno_rules[genotype]['positive_alt']) & set(mismatches))
        genotype_results[genotype]['matched_pos_kmers'] = matched
        if len(uids) > 0:
            dist =  len(mismatches) / (len(mismatches + matched))

        mutations = []
        for uid in matched:
            mutations.append(scheme_info['uid_to_mutation'][uid])
        mutations = list(set(mutations))
        ave_frac = 0
        fracs = []
        for mutation_key in mutations:
            frac = data['mutation_frac'][mutation_key]
            if frac > 0:
                fracs.append(frac)

        if len(fracs) > 0:
            ave_frac = sum(fracs) / len(fracs)

        genotype_results[genotype]['num_pos_match'] = len(matched)
        genotype_results[genotype]['num_pos_mismatch'] = len(mismatches)
        genotype_results[genotype]['kmer_genotype_dist'] = dist
        genotype_results[genotype]['ave_frac'] = ave_frac
        genotype_results[genotype]['num_pos_missing'] = missing_sites
        genotype_results[genotype]['mismatched_kmers'] = ','.join(sorted([scheme_info['uid_to_dna_name'][x] for x in mismatches]))
        genotype_results[genotype]['matched_pos_kmers_alt'] = positive_alt_match
        genotype_results[genotype]['mismatched_pos_kmers_alt'] = positive_alt_mismatch
        #del(genotype_results[genotype]['matched_pos_kmers'])

        if len(mismatches) > 0:
            genotype_results[genotype]['is_compatible'] = False

        if genotype_results[genotype]['is_compatible']:
            valid_genotypes[genotype] = genotype_results[genotype]
    df = pd.DataFrame.from_dict(genotype_results, orient='index')
    df = df.drop(columns=['matched_pos_kmers'])
    df.to_csv(outfile, header=True, sep="\t")

    return valid_genotypes

def summarize_genotype_kmers(scheme_info,kmer_results,outdir,n_threads=1):
    genotypes = scheme_info['genotypes']
    genotype_rules = scheme_info['genotype_rule_sets']

    #initialize data
    samples = list(kmer_results.keys())
    genotype_results = {}
    for sample_id in samples:
        genotype_results[sample_id] = {}
        for genotype in genotypes:
            genotype_results[sample_id][genotype] = {
                'scheme_pos_kmers_total':len(genotype_rules[genotype]['positive_uids']),
                'scheme_pos_ref_kmers_alt': len(genotype_rules[genotype]['positive_uids']),
                'kmer_genotype_dist':0,
                'num_pos_match': 0,
                'num_pos_mismatch': 0,
                'num_pos_missing': 0,
                'ave_frac':0,
                'matched_pos_kmers':[],
                'matched_pos_kmers_alt': [],
                'mismatched_pos_kmers_alt': [],
                'mismatched_kmers':[],
                'is_compatible':True
            }

    mutation_to_uids = scheme_info['mutation_to_uid']
    uid_to_mutation = scheme_info['uid_to_mutation']
    uid_to_state =  scheme_info['uid_to_state']
    num_kmers = len(uid_to_mutation)
    mutation_to_uid = {}
    for mutation_key in scheme_info['mutation_to_uid']:
        mutation_to_uid[mutation_key] = set(scheme_info['mutation_to_uid'][mutation_key])
    scheme_params = {
        'genotype_rule_sets':scheme_info['genotype_rule_sets'],
        'uid_to_mutation':scheme_info['uid_to_mutation'],
        'mutation_to_uid': mutation_to_uid,
        'uid_to_state':scheme_info['uid_to_state'],
        'uid_to_dna_name': scheme_info['uid_to_dna_name'],
        'mutation_profiles':scheme_info['mutation_profiles'],
        'mutation_positions':scheme_info['mutation_positions'],

    }

    if n_threads > 1:
        pool = Pool(processes=n_threads)

    for sample_id in samples:
        data = kmer_results[sample_id]
        outfile = os.path.join(outdir,"{}.genotypes.txt".format(sample_id))
        if n_threads == 1:
            genotype_results[sample_id] = compare_sample_to_genotypes(data,genotype_results[sample_id],scheme_params,outfile)
        else:
            genotype_results[sample_id] = pool.apply_async(compare_sample_to_genotypes,
                                                           (data,genotype_results[sample_id],scheme_params,outfile))

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

def call_compatible_genotypes(scheme_info, kmer_results, outdir, n_threads=1):
    logging.info("Comparing kmer profiles against known genotypes")
    genotype_results = summarize_genotype_kmers(scheme_info, kmer_results, outdir, n_threads)
    called_genotypes = {}
    logging.info("Identifying compatible genotypes")
    kmers_to_genotype = scheme_info['kmer_to_genotypes']
    for sample_id in genotype_results:
        valid_genotypes = {}
        parsiony_genotypes = {}
        called_genotypes[sample_id] = {'called_genotypes':genotype_results[sample_id]}
        for genotype in genotype_results[sample_id]:
            if not genotype_results[sample_id][genotype] or not genotype_results[sample_id][genotype]['is_compatible']:
                continue
            valid_genotypes[genotype] = genotype_results[sample_id][genotype]['kmer_genotype_dist']
        #order valid genotypes by lowest kmer distance
        valid_genotypes = {k: v for k, v in sorted(valid_genotypes.items(), key=lambda item: item[1])}
        genotypes_with_exclusive_kmers = {}
        pos_kmers = {}
        for genotype in valid_genotypes:
            kmers =  genotype_results[sample_id][genotype]['matched_pos_kmers']
            for k in kmers:
                if not k in pos_kmers:
                    pos_kmers[k] = []
                pos_kmers[k].append(genotype)
        for k in pos_kmers:
            if len(pos_kmers[k]) == 1:
                genotype = pos_kmers[k][0]
                if not genotype in genotypes_with_exclusive_kmers:
                    genotypes_with_exclusive_kmers[genotype] = []
                genotypes_with_exclusive_kmers[genotype].append(k)

        #Resolve kmers first by genotypes which have an exclusive kmer
        assigned_kmers = []

        for genotype in genotypes_with_exclusive_kmers:
            unassigned = list(set(genotype_results[sample_id][genotype]['matched_pos_kmers']) - set(assigned_kmers))
            if len(unassigned) == 0:
                continue
            assigned_kmers.extend(unassigned)
            parsiony_genotypes[genotype] = unassigned

        # Resolve kmers by genotypes which have the lowest distance
        prev_dist = 1
        prev_genotype = ''
        for genotype in valid_genotypes:
            dist = valid_genotypes[genotype]
            unassigned = list(set(genotype_results[sample_id][genotype]['matched_pos_kmers']) - set(assigned_kmers))
            if dist == prev_dist and prev_genotype != '':
                parsiony_genotypes[genotype] = copy.deepcopy(parsiony_genotypes[prev_genotype])
                if unassigned is not None and len(unassigned) > 0:
                    parsiony_genotypes[genotype] = parsiony_genotypes[genotype].extend(unassigned)
                continue

            if len(unassigned) == 0:
                continue

            prev_dist = dist
            prev_genotype = genotype
            assigned_kmers.extend(unassigned)
            parsiony_genotypes[genotype] = unassigned

        called_genotypes[sample_id]['called_genotypes'] = parsiony_genotypes

    return called_genotypes

def add_profile_md5(sampleManifest,scheme_info):
    #scheme lookups
    uid_to_dna_name = scheme_info['uid_to_dna_name']
    uid_to_state = scheme_info['uid_to_state']
    mutation_profiles = {}
    for sample_id in sampleManifest:
        kmer_uids = sampleManifest[sample_id]['detected_scheme_kmers']
        p = []
        for uid in kmer_uids:
            state = uid_to_state[uid]
            if state == 'ref':
                continue
            dna_name = uid_to_dna_name[uid]
            p.append(dna_name)

        md5 = calc_md5(",".join([str(x) for x in p]))
        if not md5 in mutation_profiles:
            phrase = ' '.join(generate_random_phrase())
        else:
            phrase = mutation_profiles[md5]
        sampleManifest[sample_id ]['md5'] = md5
        sampleManifest[sample_id ]['md5_phrase'] = phrase
    return sampleManifest

def write_sample_summary_results(sampleManifest,scheme_info,out_file,max_features=20):
    gene_features = scheme_info['gene_features']
    gene_features.sort()
    num_gene_features = len(gene_features)
    header = TYPER_SAMPLE_SUMMARY_HEADER_BASE
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

def write_genotype_report(sample_genotype_results,scheme_info,outfile):
    GENO_HEADER = [
        'sample_id',
        'genotype',
        'genotype_kmer_dist',
        'is_compatible',
        'num_scheme_required_kmers',
        'num_required_match_kmers',
        'num_match_alt_kmers',
        'match_alt_kmers',
        'match_alt_mutations',
        'num_required_mismatch_kmers',
        'mismatch_required_kmers',
        'mismatch_mutations'
    ]
    uid_to_state = scheme_info['uid_to_state']
    uid_to_dna_name = scheme_info['uid_to_dna_name']
    fh = open(outfile,'w')
    fh.write("{}\n".format("\t".join(GENO_HEADER)))
    for sample_id in sample_genotype_results:
        data = sample_genotype_results[sample_id]['all_genotypes']
        for genotype in data:
            alt_kmers = []
            alt_mutations = []
            for uid in data[genotype]['matched_pos_kmers']:
                if uid_to_state[uid] == 'ref':
                    continue
                alt_kmers.append(uid)
                alt_mutations.append(uid_to_dna_name[uid])

            mismatch_mutations = []
            for uid in data[genotype]['mismatched_kmers']:
                mismatch_mutations.append(uid_to_dna_name[uid])

            row = [sample_id,
                   genotype,
                   data[genotype]['kmer_genotype_dist'],
                   data[genotype]['is_compatible'],
                   data[genotype]['scheme_pos_kmers_total'],
                   data[genotype]['num_pos_match'],
                   len(data[genotype]['matched_pos_kmers_alt']),
                   "{}".format(", ".join([str(x) for x in data[genotype]['matched_pos_kmers_alt']])),
                   "{}".format(", ".join([str(x) for x in alt_mutations])),
                   len(data[genotype]['mismatched_kmers']),
                   "{}".format(", ".join([str(x) for x in data[genotype]['mismatched_kmers']])),
                   "{}".format(", ".join([str(x) for x in list(mismatch_mutations)])),
                   ]
            fh.write("{}\n".format("\t".join([str(x) for x in row])))
    fh.close()

def calc_genotype_frac(raw_kmer_counts,genotype_assigned_kmers,uid_to_mutation,mutation_to_uid,uid_to_dnaname,uid_to_state,min_cov_frac):
    genotype_fracs = {}
    for genotype in genotype_assigned_kmers:
        mutations_involved = {}
        dna_names = []
        if not isinstance(genotype_assigned_kmers[genotype],list):
            continue
        for uid in genotype_assigned_kmers[genotype]:
            mutation_key = uid_to_mutation[uid]
            freq = raw_kmer_counts[uid]
            state = uid_to_state[uid]
            if not mutation_key in mutations_involved:
                mutations_involved[mutation_key] = {'total':0,'pos': 0}
            if state == 'alt':
                mutations_involved[mutation_key]['pos']+= freq
            dna_names.append(uid_to_dnaname[uid])

        dna_names = list(set(dna_names))
        fracs = []
        for mutation_key in mutations_involved:
            for uid in mutation_to_uid[mutation_key]:
                freq = raw_kmer_counts[uid]
                mutations_involved[mutation_key]['total'] += freq

            frac = 0
            if mutations_involved[mutation_key]['total'] > 0:
                frac = mutations_involved[mutation_key]['pos']/mutations_involved[mutation_key]['total']
                if frac == 0:
                    frac = 1
                if frac >= min_cov_frac and frac <= 1 - min_cov_frac:
                    fracs.append(frac)

        ave_frac = 0
        if len(fracs) > 0:
            ave_frac = sum(fracs) / len(fracs)
        genotype_fracs[genotype] = {'ave_frac':ave_frac,'mutations':dna_names}

    return genotype_fracs

def create_sample_comparison_plots(scheme_info,kmer_results_df,labels,mds_outfile,dendro_outfile,feature_outfile):
    plots = {
        'mds':'',
        'sample_dendro':'',
        'feature_dendro':'',
        'coverage':''
    }
    #Create kmer dist matrix
    disMat = calc_sample_kmer_distMat(kmer_results_df)

    # Create MDS plot
    plots['mds'] = plot_mds(disMat, labels,mds_outfile)

    #Create Dendrogram heatmap
    d = dendrogram_visualization()
    labels = kmer_results_df.columns.tolist()
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
    for uid in alt_kmer_result_df.index.values.tolist():
        dna_names.append("{}:{}".format(uid,uid_to_dna_name[uid]))

    # Create Dendrogram heatmap of features
    plots['feature_dendro'] = d.build_dendrogram_heatmap_features(labels, dna_names, alt_kmer_result_df.T, feature_outfile)
    return plots

def create_sample_kmer_profile(kmer_results):
    return pd.DataFrame.from_dict(kmer_results,orient='index').transpose()

def calc_sample_kmer_distMat(kmer_results_df):
    return nan_compatible_kmer_pairwise_distmatrix(kmer_results_df)

def bin_scheme_targets(scheme,max_length,window_size=500):
    bins = []
    bin_mapping = {}
    for i in range(window_size,max_length+window_size,window_size):
        bin_mapping[i] = []
        bins.append(i)

    for mutation_key in scheme:
        for state in scheme[mutation_key]:
            for uid in scheme[mutation_key][state]:
                position = scheme[mutation_key][state][uid]['variant_end']
                for i in range(0, len(bins)):
                    bin = bins[i]
                    if position <= bin:
                        bin_mapping[bin].append(mutation_key)
                        bin_mapping[bin] = list(set(bin_mapping[bin] ))
                        break

    return bin_mapping

def create_coverage_plots(scheme,scheme_info,kmer_results,report_run_kmer_coverage):
    max_value = scheme_info['max_variant_positions']
    scheme_bin_mapping = bin_scheme_targets(scheme,max_value+1,window_size=int(max_value/100))
    binned_sample_data = {}
    for sample_id in kmer_results:
        binned_sample_data[sample_id] = {}
        for bin in scheme_bin_mapping:
            binned_sample_data[sample_id][bin] = []

        data = kmer_results[sample_id]['total_mutation_key_freq']
        for mutation_key in data:
            for bin in scheme_bin_mapping:
                if mutation_key in scheme_bin_mapping[bin]:
                    binned_sample_data[sample_id][bin].append(data[mutation_key])
                    break

        for bin in binned_sample_data[sample_id]:
            total = sum(binned_sample_data[sample_id][bin])
            average = 0
            if total > 0:
                average = total / len(binned_sample_data[sample_id][bin])
            binned_sample_data[sample_id][bin] = average

    binned_sample_data = pd.DataFrame.from_dict(binned_sample_data,orient='index')
    return create_heatmap(list(binned_sample_data.index.values), list(scheme_bin_mapping.keys()), binned_sample_data, report_run_kmer_coverage)

def QA_results(sampleManifest,min_coverage_depth,max_missing_sites,min_genome_size=0,max_genome_size=-1):
    for sample_id in sampleManifest:
        detected_sample_type = sampleManifest[sample_id]['detected_sample_type']
        reported_sample_type = sampleManifest[sample_id]['reported_sample_type']
        if detected_sample_type != reported_sample_type:
            sampleManifest[sample_id]['qc_messages'].append('Warning: sample type mismatch, reported:{} predicted:{}'.format(reported_sample_type,detected_sample_type))

        cov = sampleManifest[sample_id]['estimated_genome_cov']

        if cov < min_coverage_depth and sampleManifest[sample_id]['file_type'] == 'fastq':
            sampleManifest[sample_id]['qc_messages'].append(
                'Fail: low sequencing coverage {}'.format(cov))
        missing = sampleManifest[sample_id]['total_scheme_mutations'] - sampleManifest[sample_id]['detected_scheme_mutations']

        if missing > max_missing_sites:
            sampleManifest[sample_id]['qc_messages'].append(
                'Fail: low detected scheme targets {}'.format(missing))

        genomeSize = sampleManifest[sample_id]['est_genome_size']
        num_targets_found = sampleManifest[sample_id]['est_genome_size']
        if genomeSize < min_genome_size and genomeSize != 0:
            sampleManifest[sample_id]['qc_messages'].append(
                'Fail: genome size {} too small'.format(num_targets_found))
        elif genomeSize > max_genome_size and max_genome_size >= min_genome_size:
            sampleManifest[sample_id]['qc_messages'].append(
                'Fail: genome size {} too large'.format(num_targets_found))

    return sampleManifest

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
            'detected_scheme_mutations': 0,
            'detected_scheme_mixed_mutations': 0,
            'primary_genotype':'',
            'primary_genotype_frac':'',
        }

    return sampleManifest

def calc_genome_size(sampleManifest,outdir,min_cov,kLen=21,n_threads=1):
    results = {}
    if n_threads > 1:
        pool = Pool(processes=n_threads)

    for sample_id in sampleManifest:
        fileType = sampleManifest[sample_id]['file_type']
        tmpFile = os.path.join(outdir, "{}.kmc.tmp1".format(sample_id))
        tmpDir = os.path.join(outdir, "{}.kmc.tmp2".format(sample_id))
    # get GenomeSize using KMC based on fwd read or fasta
        if n_threads == 1:
            results[sample_id] = kmc_summary(sampleManifest[sample_id]['raw_seq_files'][0], tmpFile,
                               tmpDir, fileType,
                               min_cov, kLen, n_threads)
        else:
            results[sample_id] = pool.apply_async(kmc_summary,
                                                           (sampleManifest[sample_id]['raw_seq_files'][0], tmpFile,
                               tmpDir, fileType,
                               min_cov, kLen, n_threads))

    if n_threads > 1:
        pool.close()
        pool.join()

        for sample_id in results:
            results[sample_id] = results[sample_id].get()

    for sample_id in results:
        for field in results[sample_id]:
            sampleManifest[sample_id][field] = results[sample_id][field]

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

def perform_kmerSearch(sampleManifest,scheme_info,min_cov,aho,nthreads=1):
    kmer_results = {}
    # Identify kmers in each sample
    if nthreads > 1:
        pool = Pool(processes=nthreads)
    logging.info("Performing kmer searching on {} samples with {} threads".format(len(sampleManifest), nthreads))

    for sampleID in sampleManifest:
        if len(sampleManifest[sampleID]['processed_reads']) > 0:
            read_set = sampleManifest[sampleID]['processed_reads']
        else:
            read_set = sampleManifest[sampleID]['raw_seq_files']

        # Canu correct creates fasta files so need to flip this when nanopore correction has happend
        fileType = sampleManifest[sampleID]['file_type']

        if nthreads == 1:
            if fileType == 'fastq':
                kmer_results[sampleID] = perform_kmerSearch_fastq(scheme_info['uid_to_kseq'],
                                                                  scheme_info['kseq_to_uids'], aho['scheme'], read_set)
            else:
                kmer_results[sampleID] = perform_kmerSearch_fasta(scheme_info['uid_to_kseq'],
                                                                  scheme_info['kseq_to_uids'], aho['scheme'],
                                                                  read_fasta(read_set[0]), min_cov)
        else:
            if fileType == 'fastq':
                kmer_results[sampleID] = pool.apply_async(perform_kmerSearch_fastq, (
                scheme_info['uid_to_kseq'], scheme_info['kseq_to_uids'], aho['scheme'], read_set))
            else:
                kmer_results[sampleID] = pool.apply_async(perform_kmerSearch_fasta,
                                                          (scheme_info['uid_to_kseq'], scheme_info['kseq_to_uids'],
                                                           aho['scheme'], read_fasta(read_set[0]), min_cov))

    # Extract results in multithreaded mode
    if nthreads > 1:
        pool.close()
        pool.join()
        for sampleID in sampleManifest:
            kmer_results[sampleID] = kmer_results[sampleID].get()

    return kmer_results


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
    detection_limit = min_cov
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
    nthreads = int(cmd_args.n_threads)
    strain_profiles = cmd_args.strain_profiles
    primers = cmd_args.primers
    no_template_control = cmd_args.no_template_control
    positive_control = cmd_args.positive_control
    genotype_dist_cutoff = cmd_args.genotype_dist_cutoff
    min_genome_size = cmd_args.min_genome_size
    max_genome_size = cmd_args.max_genome_size
    min_genome_cov_depth = cmd_args.min_genome_cov_depth
    max_features = cmd_args.max_features
    type_only = cmd_args.typer_only

    #Initialize scheme
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

    # output filenames
    report_run_info_log = open(os.path.join(outdir, "{}.run.info.log".format(prefix)), 'w')
    report_sample_composition_summary = os.path.join(outdir, "{}.sample_composition.report.summary.txt".format(prefix))
    report_genotype_targets = os.path.join(outdir, "{}.sample.genotype.targets.txt".format(prefix))
    report_sample_kmer_profiles = os.path.join(outdir, "{}.sample.kmer.profiles.txt".format(prefix))
    report_run_metrics = open(os.path.join(outdir, "{}.run.metrics".format(prefix)), 'w')
    report_summary = os.path.join(outdir, "{}.run.summary.html".format(prefix))


    #Plots
    report_run_kmer_dendrogram = os.path.join(outdir, "{}.run_kmer_dendrogram.html".format(prefix))
    report_run_kmer_heatmap = os.path.join(outdir, "{}.run_kmer_heatmap.html".format(prefix))
    report_run_kmer_mds = os.path.join(outdir, "{}.run_kmer_mds.html".format(prefix))
    report_run_kmer_coverage = os.path.join(outdir, "{}.run_kmer_coverage.html".format(prefix))

    #Logging
    analysis_date = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    report_run_info_log.write("Start Time\t{}\n".format(analysis_date))
    report_run_info_log.write("Output Directory\t{}\n".format(outdir))
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
    logger.info("Determine read stats for {} samples".format(len(seqManifest)))

    fastp_dir = os.path.join(outdir, "fastp")
    if not os.path.isdir(fastp_dir) and not type_only:
        logger.info("Creating processed read directory {}".format(fastp_dir))
        os.mkdir(fastp_dir, 0o755)
    else:
        logger.info(
            "Results directory {} already exits, will overwrite any results files here".format(fastp_dir))

    if len(seqManifest) == 0:
        logger.error("Error no fasta/fastq files to process")
        sys.exit()

    sampleManifest = init_sample_manifest(seqManifest, scheme_info, scheme_name, analysis_date, sample_type, seqTech)

    #process reads
    if not type_only:
        sampleManifest = calc_genome_size(sampleManifest, outdir, min_cov, kLen=21, n_threads=nthreads)
        sampleManifest = process_reads(sampleManifest,perform_read_correction,read_dir,seqTech,fastp_dir,min_read_len,trim_front_bp,trim_tail_bp,perform_read_dedup,n_threads=nthreads)

    kmer_results = perform_kmerSearch(sampleManifest,scheme_info,min_cov,aho,nthreads)

    #identify contamination from no template control
    contaminants = {}
    if no_template_control is not None and no_template_control in seqManifest:
        logger.info("Processing negative control: {}".format(no_template_control))
        fileType = seqManifest[no_template_control]['seq_type']
        read_set = seqManifest[no_template_control]['seq_files']
        if fileType == 'fastq':
            contaminants = perform_kmerSearch_fastq(scheme_info['uid_to_kseq'], scheme_info['kseq_to_uids'],
                                                              aho['scheme'], read_set)
        else:
            contaminants = perform_kmerSearch_fasta(scheme_info['uid_to_kseq'], scheme_info['kseq_to_uids'],
                                                              aho['scheme'], read_fasta(read_set[0]), min_cov)
    found_no_template_kmers = {}
    for uid in contaminants:
        freq = contaminants[uid]
        if freq >= min_cov:
            found_no_template_kmers[uid] = freq

    found_no_template_uids = set(list(found_no_template_kmers.keys()))
    logger.info("Found {} kmers in no template control".format(len(found_no_template_uids)))
    if len(found_no_template_uids) > 0:
        logger.info("Performing no template subtraction")

    missing_mutations = {}
    #subtract contaminants
    for sampleID in kmer_results:
        if sampleID == no_template_control:
            continue
        temp_freq = {}
        for mutation_key in scheme_info['mutation_to_uid']:
            mutation_uids = scheme_info['mutation_to_uid'][mutation_key]
            found_uids = []
            freqs = []
            alt_freqs = []
            for uid in mutation_uids:
                if uid in kmer_results[sampleID]:
                    freq =kmer_results[sampleID][uid]
                    freqs.append(freq)
                    if freq >= min_cov:
                        found_uids.append(uid)
                        if scheme_info['uid_to_state'] == 'alt':
                            alt_freqs.append(freq)
            if len(found_uids) == 0:
                if not mutation_key in missing_mutations:
                    missing_mutations[mutation_key] = 0
                missing_mutations[mutation_key] += 1
            noContam = set(found_uids) - found_no_template_uids
            if len(noContam) == 0:
                continue
            ovContam = list(set(found_uids) & found_no_template_uids)
            for uid in ovContam:
                conFreq = contaminants[uid]
                sampleFreq = kmer_results[sampleID][uid]
                newFreq = sampleFreq - conFreq
                if newFreq < 0:
                    newFreq = 0
                kmer_results[sampleID][uid] = newFreq

    logger.info("Kmer searching complete")
    kmer_results_df = create_sample_kmer_profile(kmer_results)
    kmer_results_df.to_csv(report_sample_kmer_profiles,sep="\t",header=True)


    plots = {
        'mds':None,
        'sample_dendro':None,
        'feature_dendro':None,
        'coverage':None
    }

    # create a plot of sample similarity for a multi-sample run
    if len(sampleManifest)  > 1 and no_plots == False:
        if len(sampleManifest) <= 1000 or force_plots:
            logging.info('Creating sample comparison plots')
            plots = create_sample_comparison_plots(scheme_info,kmer_results_df, kmer_results_df.columns.tolist(), report_run_kmer_mds, report_run_kmer_dendrogram,report_run_kmer_heatmap)


    #process kmer results
    logger.info("Processing kmer results")
    kmer_results = process_kmer_results(scheme_info, kmer_results, min_cov, min_cov_frac)

    logger.info("Creating coverage plot")
    plots['coverage'] = create_coverage_plots(scheme, scheme_info, kmer_results, report_run_kmer_coverage)

    mixed_muations = {}
    for sample_id in kmer_results:
        sampleManifest[sample_id]['detected_scheme_mixed_mutations'] = kmer_results[sample_id]['num_mixed_sites']
        sampleManifest[sample_id]['ave_scheme_kmers_freq'] = kmer_results[sample_id]['average_kmer_freq']
        sampleManifest[sample_id]['detected_scheme_mutations'] = sampleManifest[sample_id]['total_scheme_mutations'] - kmer_results[sample_id]['num_missing_sites']
        sampleManifest[sample_id]['detected_scheme_kmers'] = kmer_results[sample_id]['detected_scheme_kmers']
        sampleManifest[sample_id]['num_detected_scheme_kmers'] = len(kmer_results[sample_id]['detected_scheme_kmers'])
        if sampleManifest[sample_id]['detected_scheme_mixed_mutations']  > max_mixed_sites :
            sampleManifest[sample_id]['detected_sample_type'] = 'multi'
        for mutation_key in kmer_results[sample_id]['mixed_sites']:
            if not mutation_key in mixed_muations:
                missing_mutations[mutation_key] = 0
            missing_mutations[mutation_key] += 1

    #Add MD5 calculation and MD5 phrase for duplicate sample identification
    sampleManifest = add_profile_md5(sampleManifest, scheme_info)

    #Compare kmer results with genotype rules in the scheme
    logger.info("Comparing kmer results to genotypes")
    genotype_dir = os.path.join(outdir,"genotype.reports")
    if not os.path.isdir(genotype_dir):
        logger.info("Creating processed read directory {}".format(genotype_dir))
        os.mkdir(genotype_dir, 0o755)
    else:
        logger.info("Results directory {} already exits, will overwrite any results files here".format(genotype_dir))

    sample_genotype_results = call_compatible_genotypes(scheme_info, kmer_results, genotype_dir,  nthreads)

    #Assign primary genotype if possible
    for sample_id in sample_genotype_results:
        if len(sample_genotype_results[sample_id]['called_genotypes']) == 1 :
            sampleManifest[sample_id]['primary_genotype'] = list(sample_genotype_results[sample_id]['called_genotypes'].keys())[0]
            sampleManifest[sample_id]['primary_genotype_frac'] = 1.0
        else:
            genotype_mutation_fracs = calc_genotype_frac(kmer_results[sample_id]['raw_kmer_freq'],sample_genotype_results[sample_id]['called_genotypes'],
                                                         scheme_info['uid_to_mutation'],
                                                         scheme_info['mutation_to_uid'],
                                                         scheme_info['uid_to_dna_name'],
                                                         scheme_info['uid_to_state'],min_cov_frac)
            pGenotype = []
            pGenotype_frac = - 1
            for geno in genotype_mutation_fracs:
                frac = genotype_mutation_fracs[geno]['ave_frac']
                if frac > pGenotype_frac:
                    pGenotype_frac = frac
                    pGenotype = [geno]
                elif frac == pGenotype_frac:
                    pGenotype.append(pGenotype)
            if len(pGenotype) != 1:
                pGenotype = 'n/a'
                pGenotype_frac = 0
            else:
                pGenotype = pGenotype[0]

            sampleManifest[sample_id]['primary_genotype'] = pGenotype
            sampleManifest[sample_id]['primary_genotype_frac'] = pGenotype_frac
        sampleManifest[sample_id]['compatible_genotypes'] = list(sample_genotype_results[sample_id]['called_genotypes'].keys())

    sampleManifest = QA_results(sampleManifest, min_genome_cov_depth, max_missing_sites, min_genome_size, max_genome_size)

    sample_type_mismatch = 0
    num_fail_samples = 0
    ave_kmer_freq_list = []
    genotypes_found = {}
    failed_sample_ids = []

    for sample_id in sample_genotype_results:
        sType = sampleManifest[sample_id]['reported_sample_type']
        pType = sampleManifest[sample_id]['detected_sample_type']
        if sType != pType:
            sample_type_mismatch+=1
        if len(sampleManifest[sample_id]['qc_messages']) > 0:
            num_fail_samples+=1
            failed_sample_ids.append(sample_id)
        ave_kmer_freq_list.append(sampleManifest[sample_id]['ave_scheme_kmers_freq'])
        compatible_genotypes = sampleManifest[sample_id]['compatible_genotypes']
        for g in compatible_genotypes:
            if not g in genotypes_found:
                genotypes_found[g] =0
            genotypes_found[g] += 1

    kmer_freq_ave = 0
    kmer_freq_stdev = 0
    if len(ave_kmer_freq_list) >= 1:
        kmer_freq_ave= sum(ave_kmer_freq_list) / len(ave_kmer_freq_list)
        if len(ave_kmer_freq_list) > 1:
            kmer_freq_stdev = statistics.stdev(ave_kmer_freq_list)


    #create genotype plot
    logger.info("Creating genotype abundance chart")
    num_genotypes_found = len(genotypes_found)
    genotypes_found = {k: v for k, v in sorted(genotypes_found.items(), key=lambda item: item[1], reverse=True)}
    genotypes_found_df = pd.DataFrame.from_dict(genotypes_found, orient='index', columns=['count'])
    genotypes_found_df['genotype'] = list(genotypes_found_df.index)
    fig = px.bar(genotypes_found_df, x='genotype', y='count')
    fig.update_layout({'width': 1000,
                       'height': 800,
                       'showlegend': False,
                       'hovermode': 'closest',
                       })
    plots['genotype_abundance'] = fig

    # create mutation plot
    logger.info("Creating missing target chart")
    missing_mutations = {k: v for k, v in sorted(missing_mutations.items(), key=lambda item: item[1], reverse=True)}
    missing_mutations_df = pd.DataFrame.from_dict(missing_mutations, orient='index', columns=['count'])
    missing_mutations_df['mutation'] = list(missing_mutations_df.index)
    fig = px.bar(missing_mutations_df, x='mutation', y='count')
    fig.update_layout({'width': 1000,
                       'height': 800,
                       'showlegend': False,
                       'hovermode': 'closest',
                       })
    plots['missing_features'] = fig

    # create mutation plot
    logger.info("Creating mixed target summary chart")
    mixed_muations = {k: v for k, v in sorted(mixed_muations.items(), key=lambda item: item[1], reverse=True)}
    mixed_muations_df = pd.DataFrame.from_dict(mixed_muations, orient='index', columns=['count'])
    mixed_muations_df['mutation'] = list(mixed_muations_df.index)
    fig = px.bar(mixed_muations_df, x='mutation', y='count')
    fig.update_layout({'width': 1000,
                       'height': 800,
                       'showlegend': False,
                       'hovermode': 'closest',
                       })
    plots['coverage_mixed'] = fig

    # create sample summary table
    write_sample_summary_results(sampleManifest, scheme_info, report_sample_composition_summary, max_features)

    #create batch html report if more than one sample present
    if len(sampleManifest) > 1 and not no_plots:
        template_html = Path(BATCH_HTML_REPORT).read_text()
        html_report_data = {
            'analysis_date': analysis_date,
            'total_samples': len(sampleManifest),
            'positive_control_id':positive_control,
            'negative_control_id': no_template_control,
            'sample_type': type,
            'sample_type_mismatch':sample_type_mismatch,
            'strain_scheme_name':scheme_name,
            'num_targets_detected':len(scheme_info['mutation_to_uid']) - len(missing_mutations),
            'ave_kmer_freq':  kmer_freq_ave,
            'stdev_kmer_freq': kmer_freq_stdev,
            'count_failed_samples': num_fail_samples,
            'failed_sample_ids': ', '.join([str(x) for x in failed_sample_ids]),

            'coverage_plot': plot(figure_or_data=plots['coverage'],include_plotlyjs=False,
                 output_type='div'),
            'coverage_plot_caption': FIGURE_CAPTIONS['coverage_plot_caption'],

            'mixed_target_plot': plot(figure_or_data=plots['coverage_mixed'],include_plotlyjs=False,
                 output_type='div'),

            'mixed_target_plot_caption': FIGURE_CAPTIONS['mixed_target_plot_caption'],

            'missing_target_plot': plot(figure_or_data=plots['missing_features'],include_plotlyjs=False,
                 output_type='div'),

            'missing_target_plot_caption': FIGURE_CAPTIONS['missing_features_plot_caption'],
            'positive_control_found_targets': 0,
            'positive_control_missing_targets':0,
            'positive_control_mixed_targets': 0,
            'negative_control_found_targets': len(found_no_template_uids),
            'negative_control_missing_targets':0,
            'negative_control_mixed_targets': 0,

            'num_genotypes_found': num_genotypes_found,
            'genotypes_found': ', '.join([str(x) for x in genotypes_found]),
            'genotype_abundance_plot': plot(figure_or_data=plots['genotype_abundance'],include_plotlyjs=False,
                 output_type='div'),
            'genotype_abundance_plot_caption': FIGURE_CAPTIONS['genotype_abundance_plot_caption'],

            'sample_mds':plot(figure_or_data=plots['mds'],include_plotlyjs=False,
                 output_type='div'),
            'sample_mds_caption':FIGURE_CAPTIONS['sample_mds_caption'],

            'sample_dendro':plot(figure_or_data=plots['sample_dendro'],include_plotlyjs=False,
                 output_type='div'),
            'sample_dendro_caption': FIGURE_CAPTIONS['sample_dendro_caption'],

            'feature_dendro':plot(figure_or_data=plots['feature_dendro'],include_plotlyjs=False,
                 output_type='div'),
            'feature_dendro_caption': FIGURE_CAPTIONS['feature_dendro_caption'],

            'analysis_parameters': analysis_parameters

        }
        batch_sample_html_report(template_html, {'data':html_report_data}, report_summary)
    elif len(sampleManifest) > 1 :
        template_html = Path(BATCH_HTML_REPORT).read_text()
        html_report_data = {
            'analysis_date': analysis_date,
            'total_samples': len(sampleManifest),
            'positive_control_id': positive_control,
            'negative_control_id': no_template_control,
            'sample_type': type,
            'sample_type_mismatch': sample_type_mismatch,
            'strain_scheme_name': scheme_name,
            'num_targets_detected': len(scheme_info['mutation_to_uid']) - len(missing_mutations),
            'ave_kmer_freq': kmer_freq_ave,
            'stdev_kmer_freq': kmer_freq_stdev,
            'count_failed_samples': num_fail_samples,
            'failed_sample_ids': ', '.join([str(x) for x in failed_sample_ids]),

            'coverage_plot': plot(figure_or_data=plots['coverage'], include_plotlyjs=False,
                                  output_type='div'),
            'coverage_plot_caption': FIGURE_CAPTIONS['coverage_plot_caption'],

            'mixed_target_plot': plot(figure_or_data=plots['coverage_mixed'], include_plotlyjs=False,
                                      output_type='div'),

            'mixed_target_plot_caption': FIGURE_CAPTIONS['mixed_target_plot_caption'],

            'missing_target_plot': plot(figure_or_data=plots['missing_features'], include_plotlyjs=False,
                                        output_type='div'),

            'missing_target_plot_caption': FIGURE_CAPTIONS['missing_features_plot_caption'],
            'positive_control_found_targets': 0,
            'positive_control_missing_targets': 0,
            'positive_control_mixed_targets': 0,
            'negative_control_found_targets': len(found_no_template_uids),
            'negative_control_missing_targets': 0,
            'negative_control_mixed_targets': 0,

            'num_genotypes_found': num_genotypes_found,
            'genotypes_found': ', '.join([str(x) for x in genotypes_found]),

            'analysis_parameters': analysis_parameters

        }
        batch_sample_html_report(template_html, {'data': html_report_data}, report_summary)

