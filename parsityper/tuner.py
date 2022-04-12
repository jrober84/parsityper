#!/usr/bin/python
import logging
import multiprocessing as mp
import os, time
import sys
from argparse import (ArgumentParser)

import pandas as pd
from scipy.stats import entropy
from parsityper.helpers import init_console_logger, read_tsv
from parsityper.scheme import parseScheme, constructSchemeLookups
from parsityper.version import __version__

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Parsityper genotype tuning')
    parser.add_argument('--input_meta', type=str, required=True,
                        help='TSV file of sample_id,genotype')
    parser.add_argument('--input_profile', type=str, nargs='+',required=True,
                        help='Pre-computed kmer result profile for samples')
    parser.add_argument('--scheme', type=str, required=True,
                        help='TSV formated kmer scheme')
    parser.add_argument('--outdir', type=str, required=True,
                        help='output directory')
    parser.add_argument('--prefix', type=str, required=False,
                        help='output file prefix',default='parsityper')
    parser.add_argument('--update', required=False,
                        help='Do not blank existing rules, only update rules based on supplied samples',action='store_true')
    parser.add_argument('--no_mixed', required=False,
                        help='Remove sites where both ref and alt exist in the same sample',action='store_true')
    parser.add_argument('--min_cov', type=int, required=False,
                        help='Absolute minimum kmer coverage for fastq read detection default=auto determine coverage',default=50)
    parser.add_argument('--min_cov_frac', type=float, required=False,
                        help='Minimum percentage of total pool required for k-mer detection range 0 - 1.0 (default=0.05)', default=0.05)
    parser.add_argument('--min_alt_frac', type=float, required=False,
                        help='Minimum fraction of isolates positive for mutation for it to be positive 0 - 1.0 (default=0.1)', default=0.95)
    parser.add_argument('--min_ref_frac', type=float, required=False,
                        help='Minimum fraction of isolates positive for mutation for it to be positive 0 - 1.0 (default=0.1)', default=1)
    parser.add_argument('--min_alt_freq', type=int, required=False,
                        help='Minimum number of isolates positive for mutation for it to be valid for the scheme (default=1)', default=0)
    parser.add_argument('--max_frac_missing', type=float, required=False,
                        help='Maximum number of sequences allowed to be missing a mutation for it to stay included', default=0.25)
    parser.add_argument('--max_conflict', type=float, required=False,
                        help='Maximum percentage of genotype samples before droping as a requirement range 0 - 1.0 (default=0.05)', default=0.05)
    parser.add_argument('--n_threads', type=int, required=False,
                        help='output directory',default=1)
    parser.add_argument('-V', '--version', action='version', version='%(prog)s {}'.format(__version__))

    return parser.parse_args()

def get_genotype_mapping(metadata_df):
    '''
    :param metadata_df:
    :type metadata_df:
    :return:
    :rtype:
    '''
    mapping = {}
    for row in metadata_df.itertuples():
        sample_id = row.sample_id
        genotype = row.genotype
        mapping[sample_id] = genotype
    return mapping

def identify_mixed_sites(sample_mapping,kmer_file,mixed_mutations,scheme_info,min_cov):
    select_uids = []
    select_geno = []
    for uid in mixed_mutations:
        select_uids.append(uid)
        select_geno += mixed_mutations[uid]
    select_uids = sorted(list(set(select_uids)))
    select_geno = sorted(list(set(select_uids)))


    kmer_FH = open(kmer_file, 'r')
    samples = next(kmer_FH).rstrip().split("\t")
    valid_cols_range = []
    sample_kcounts = {}
    for i in range(0, len(samples)):
        sample_id = samples[i]
        if sample_id in sample_mapping:
            genotype = sample_mapping[sample_id]
            if genotype in select_geno:
                valid_cols_range.append(i)
                sample_kcounts[sample_id] = {}
                for uid in select_uids:
                    sample_kcounts[sample_id][uid] = 0
    uid = 0
    for line in kmer_FH:
        #skip rows not involved in potential conflicts
        if uid not in select_uids:
            uid+=1
            continue
        row = line.rstrip().split("\t")
        for i in valid_cols_range:
            sample_id = samples[i]
            value = int(row[i])
            if value >= min_cov:
                sample_kcounts[sample_id][uid] = value
        uid += 1
    kmer_FH.close()

    stime = time.time()
    print("Calculating frac {}".format(time.time() - stime))
    sample_mixed_sites = set()
    for sample_id in sample_kcounts:
        for uid in mixed_mutations:
            if uid in sample_mixed_sites:
                continue
            mutation = scheme_info['uid_to_mutation'][uid]
            uids = scheme_info['mutation_to_uid'][mutation]
            is_alt_present = False
            is_ref_present = False
            for uid in uids:
                value = sample_kcounts[sample_id][uid]
                if value < min_cov:
                    continue
                if scheme_info['uid_to_state'][uid] == 'alt':
                    is_alt_present = True
                else:
                    is_ref_present = True
            if is_ref_present and is_alt_present:
                sample_mixed_sites = sample_mixed_sites | set(uids)
    return sample_mixed_sites

def get_genotype_kmer_counts(sample_mapping,kmer_file,scheme_info,num_kmers,min_cov=1):
    genotypes = sorted(list(set(list(sample_mapping.values()))))
    genotype_kmer_counts = {}
    for genotype in genotypes:
        genotype_kmer_counts[genotype] = [0] * num_kmers

    kmer_counts = {}
    for i in range(0, num_kmers):
        kmer_counts[i] = 0
    kmer_FH = open(kmer_file, 'r')
    samples = next(kmer_FH).rstrip().split("\t")
    valid_cols_range = []
    genotype_cols = [0] * len(samples)
    for i in range(0, len(samples)):
        sample_id = samples[i]
        if sample_id in sample_mapping:
            valid_cols_range.append(i)
            genotype_cols[i] = sample_mapping[sample_id]
    stime = time.time()
    uid = 0
    for line in kmer_FH:
        row = line.rstrip().split("\t")
        for i in valid_cols_range:
            genotype = genotype_cols[i]
            value = int(row[i])
            if value >= min_cov:
                genotype_kmer_counts[genotype][uid] += 1
                kmer_counts[uid] += 1
        uid += 1
    kmer_FH.close()

    if num_kmers != uid:
        logging.error(
            "The number of rows in the profile do not match the number of scheme kmers. Check that it is not truncated")
        logging.error(
            "File: {}, Scheme Kmers:{}, Profile Kmers:{}".format(kmer_file, num_kmers, uid))
        sys.exit()

    #Look for potentially mixed sites based on genotype counts
    geno_fracs = calc_site_frac(genotype_kmer_counts, scheme_info)
    mutation_range = range(0, len(scheme_info['mutation_to_uid']))
    mixed_mutations = {}
    for genotype in geno_fracs:
        for mutation in mutation_range:
            value = geno_fracs[genotype][mutation]
            if value > 0 and value < 1:
                if not mutation in mixed_mutations:
                    mixed_mutations[mutation] = []
                mixed_mutations[mutation].append(genotype)

    mixed_uids = set()
    if len(mixed_mutations) > 0:
        mixed_uids = identify_mixed_sites(sample_mapping,kmer_file,mixed_mutations,scheme_info,min_cov)

    return {'genotype_kmer_counts': genotype_kmer_counts, 'mixed_sites': mixed_uids}


def get_genotype_kmer_counts_bck(sample_mapping,kmer_file,scheme_info,num_kmers,min_cov=1):
    genotypes = sorted(list(set(list(sample_mapping.values()))))
    genotype_kmer_counts = {}
    for genotype in genotypes:
        genotype_kmer_counts[genotype] = [0] * num_kmers
    sample_kcounts = {}
    for sample_id in sample_mapping:
        sample_kcounts[sample_id] = [0] * num_kmers

    kmer_counts = {}
    for i in range(0,num_kmers):
        kmer_counts[i] = 0
    kmer_FH = open(kmer_file,'r')
    samples = next(kmer_FH).rstrip().split("\t")
    valid_cols_range = []
    genotype_cols = []
    for i in range(0,len(samples)):
        sample_id = samples[i]
        if sample_id in sample_mapping:
            valid_cols_range.append(i)
            genoytpe_cols.append(sample_mapping[sample_id])
    stime = time.time()
    uid = 0
    for line in kmer_FH:
        row = line.rstrip().split("\t")
        for i in valid_cols_range:
            sample_id = samples[i]
            genotype = sample_mapping[sample_id]
            value = int(row[i])
            if value >= min_cov:
                sample_kcounts[sample_id][uid] = value
                genotype_kmer_counts[genotype][uid]+=1
                kmer_counts[uid]+=1
        uid+=1
    kmer_FH.close()

    if num_kmers != uid:
        logging.error("The number of rows in the profile do not match the number of scheme kmers. Check that it is not truncated")
        logging.error(
            "File: {}, Scheme Kmers:{}, Profile Kmers:{}".format(kmer_file,num_kmers,uid))
        sys.exit()
    stime = time.time()
    fracs = calc_site_frac(sample_kcounts, scheme_info)


    mutation_range = range(0,len(scheme_info['mutation_to_uid']))
    mixed_sites = {}
    for sample_id in fracs:
        for uid in mutation_range:
            value = fracs[sample_id][uid]
            if value > 0 and value < 1:
                if not uid in mixed_sites:
                    mixed_sites[uid] = 0
                mixed_sites[uid]+=1

    return {'genotype_kmer_counts':genotype_kmer_counts,'mixed_sites':mixed_sites}

def get_genotype_counts(sample_mapping):
    counts = {}
    genotypes = sorted(list(set(list(sample_mapping.values()))))
    for genotype in genotypes:
        counts[genotype] = 0

    for sample_id in sample_mapping:
        genotype = sample_mapping[sample_id]
        counts[genotype]+=1

    return counts

def calc_shanon_entropy(value_list):
    total = sum(value_list)
    if total == 0:
        return -1
    values = []
    for v in value_list:
        values.append(v / total)
    return entropy(values)

def determine_genotype_kmer_assoc(genotype_kmer_counts,genotype_counts,scheme_info,min_ref_frac,min_alt_frac,min_alt_freq):
    rules = {}
    sEntropy = {}
    min_par_frac = max([1 - min_alt_frac, 1 - min_ref_frac])
    uid_range = range(0,len(scheme_info['uid_to_state']))
    uid_geno_counts = {}
    for uid in uid_range:
        uid_geno_counts[uid] = {}
        for genotype in genotype_kmer_counts:
            uid_geno_counts[uid][genotype] = 0
    stime = time.time()
    for genotype in genotype_kmer_counts:
        if not genotype in genotype_counts:
            continue
        genotype_count = genotype_counts[genotype]
        rules[genotype] = {'positive_uids':[],
                           'positive_ref':[],
                           'positive_alt':[],
                           'partial_uids':[],
                           'partial_ref':[],
                           'partial_alt':[]}
        for uid in uid_range:
            value = genotype_kmer_counts[genotype][uid]
            if value == 0:
                continue
            uid_geno_counts[uid][genotype] = value
            state = scheme_info['uid_to_state'][uid]

            frac = value / genotype_count

            if frac < min_par_frac:
                continue
            if state == 'ref' and frac >= min_ref_frac:
                rules[genotype]['positive_uids'].append(uid)
                rules[genotype]['positive_ref'].append(uid)
            elif state == 'alt' and frac >= min_alt_frac:
                rules[genotype]['positive_uids'].append(uid)
                rules[genotype]['positive_alt'].append(uid)
            elif state == 'ref' and frac >= min_par_frac:
                rules[genotype]['partial_uids'].append(uid)
                rules[genotype]['positive_ref'].append(uid)
            elif state == 'alt' and frac >= min_par_frac:
                rules[genotype]['partial_uids'].append(uid)
                rules[genotype]['partial_alt'].append(uid)
    print(time.time() - stime)
    kmer_rules = {}
    for uid in scheme_info['uid_to_state']:
        kmer_rules[uid] = {
            'positive_genotypes': [],
            'partial_genotypes':[]
        }
    stime = time.time()
    for genotype in rules:
        for uid in rules[genotype]['positive_uids']:
            kmer_rules[uid]['positive_genotypes'].append(genotype)
        for uid in rules[genotype]['partial_uids']:
            kmer_rules[uid]['partial_genotypes'].append(genotype)
    print(time.time() - stime)
    stime = time.time()
    for uid in uid_geno_counts:
        sEntropy[uid] = calc_shanon_entropy( list(uid_geno_counts[uid].values()))
    print(time.time() - stime)
    return {'geno_rules':rules,'entropy':sEntropy,'kmer_rules':kmer_rules}

def blank_invalid_rules(kmer_rules,num_genotypes,scheme_info):
    for mutation_key in scheme_info['mutation_to_uid']:
        uids = scheme_info['mutation_to_uid'][mutation_key]
        ref_geno = []
        alt_geno = []
        for uid in uids:
            state = scheme_info['uid_to_state'][uid]
            if state == 'ref':
                ref_geno += kmer_rules[uid]['positive_genotypes'] + kmer_rules[uid]['partial_genotypes']
                continue
            alt_geno += kmer_rules[uid]['positive_genotypes'] + kmer_rules[uid]['partial_genotypes']
            count_alt = 0
            if len(kmer_rules[uid]['positive_genotypes']) == 0:
                count_alt += 1

        if count_alt == 0 or len(set(ref_geno)) == num_genotypes:
            for uid in uids:
                kmer_rules[uid]['positive_genotypes'] = []
                kmer_rules[uid]['partial_genotypes'] = []

    return kmer_rules

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


def filter_missing_sites(kmer_counts,scheme_info,max_missing_count):
    invalid_kmers = []
    for mutation_key in scheme_info['mutation_to_uid']:
        uids = scheme_info['mutation_to_uid'][mutation_key]
        count_present = 0
        for uid in uids:
             count_present+= kmer_counts[uid]
        if count_present < max_missing_count:
            invalid_kmers += uids
    return invalid_kmers

def update_scheme(input_scheme,output_scheme,valid_uids,kmer_geno_rules,kmer_entropies):
    df = read_tsv(input_scheme)
    positive_genotypes = []
    partial_genotypes = []
    entropies = []
    for uid in kmer_geno_rules:
        positive_genotypes.append(','.join(kmer_geno_rules[uid]['positive_genotypes']))
        partial_genotypes.append(','.join(kmer_geno_rules[uid]['partial_genotypes']))
        if uid in kmer_entropies:
            entropies.append(kmer_entropies[uid])
        else:
            entropies.append(-1)

    df['positive_genotypes'] = positive_genotypes
    df['partial_genotypes'] = partial_genotypes
    df['entropy'] = entropies
    df = df[df['key'].isin(valid_uids)].reset_index(drop=True)
    df['key'] = df.index
    df.to_csv(output_scheme,sep="\t",header=True,index=False)

def run():
    # input parameters
    cmd_args = parse_args()
    input_meta = cmd_args.input_meta
    input_profile = cmd_args.input_profile
    scheme_file = cmd_args.scheme
    prefix = cmd_args.prefix
    min_cov = int(cmd_args.min_cov)
    only_update = cmd_args.update
    min_alt_frac = cmd_args.min_alt_frac
    min_ref_frac = cmd_args.min_ref_frac
    min_alt_freq = cmd_args.min_alt_freq
    no_mixed = cmd_args.no_mixed

    max_frac_missing = cmd_args.max_frac_missing
    outdir = cmd_args.outdir
    logger = init_console_logger(2)


    #result files
    scheme_outfile = os.path.join(outdir,"{}-scheme.txt".format(prefix))

    # initialize analysis directory
    if not os.path.isdir(outdir):
        logger.info("Creating analysis results directory {}".format(outdir))
        os.mkdir(outdir, 0o755)
    else:
        logger.info("Results directory {} already exits, will overwrite any results files here".format(outdir))

    # read sample metadata
    if not os.path.isfile(input_meta):
        logger.error("Error metadata file {} does not exist".format(input_meta))
        sys.exit()
    logging.info("Reading genotype assignments from {}".format(input_meta))
    metadata_df = read_tsv(input_meta)
    logging.info("Found {} lines in {}".format(len(metadata_df), input_meta))
    metadata_df['sample_id'] = metadata_df['sample_id'].astype(str)
    metadata_df['genotype'] = metadata_df['genotype'].astype(str)
    sample_mapping = get_genotype_mapping(metadata_df)


    #Read scheme
    logger.info("Initializing scheme data structure from {}".format(scheme_file))
    if not os.path.isfile(scheme_file):
        logger.error("Error scheme file {} does not exist".format(scheme_file))
        sys.exit()
    scheme = parseScheme(scheme_file)
    scheme_info = constructSchemeLookups(scheme)

    # Read Kmer header
    profile_samples = []
    for file in input_profile:
        if not os.path.isfile(file):
            logger.error("Error profile {} does not exist".format(file))
            sys.exit()

        kmer_FH = open(file,'r')
        header = next(kmer_FH).rstrip().split("\t")
        kmer_FH.close()
        profile_samples += header[1:]

    filt = {}
    for sample_id in profile_samples:
        if sample_id in sample_mapping:
            filt[sample_id] = sample_mapping[sample_id]
    sample_mapping = filt
    genotype_counts = get_genotype_counts(sample_mapping)
    num_genotypes = len(genotype_counts)
    max_missing_count = int(max_frac_missing * len(sample_mapping))

    #Confirm all samples in profile exist in metadata
    ovl = set(profile_samples) & set(list(sample_mapping.keys()))
    if len(ovl) != len(set(profile_samples)):
        logger.error("Sample id's in profile are not present in metadata {}, these will be exluded from analysis".format(set(profile_samples) - set(list(sample_mapping.keys()))))

    num_kmers = len(scheme_info['uid_to_state'])
    logger.info("Found {} samples in all profiles".format(len(profile_samples)))
    kdata = {'genotype_kmer_counts':{},'kmer_counts':{},'mixed_sites':set()}
    for genotype in genotype_counts:
        kdata['genotype_kmer_counts'][genotype] = [0] * num_kmers
    for uid in range(0,num_kmers):
        kdata['kmer_counts'][uid] = 0


    stime = time.time()
    for file in input_profile:
        logging.info("Reading profile from {}".format(file))
        sys.stdout.flush()
        data = get_genotype_kmer_counts(sample_mapping, file, scheme_info, num_kmers, min_cov)
        for genotype in data['genotype_kmer_counts']:
            for uid,value in enumerate(data['genotype_kmer_counts'][genotype]):
                kdata['genotype_kmer_counts'][genotype][uid]+= value
                kdata['kmer_counts'][uid]+=value
        kdata['mixed_sites'] = kdata['mixed_sites'] | data['mixed_sites']

    out_str = []
    for uid in kdata['mixed_sites']:
        row = []
        row.append(str(uid))
        row.append(str(scheme_info['uid_to_mutation'][uid]))
        row.append(str(scheme_info['uid_to_dna_name'][uid]))
        row.append(str(scheme_info['uid_to_aa_name'][uid]))
        row.append(str(scheme_info['uid_to_gene_feature'][uid]))
        row.append(str(scheme_info['uid_to_state'][uid]))
        row.append(str(scheme_info['uid_to_kseq'][uid]))
        out_str.append("{}".format("\t".join(row)))
    out_fh = open(os.path.join(outdir, "{}-mixed-kmers.txt".format(prefix)), 'w')
    out_fh.write("{}".format("\n".join(out_str)))
    out_fh.close()


    logger.info("Associating kmers")
    assoc_data = determine_genotype_kmer_assoc(kdata['genotype_kmer_counts'], genotype_counts, scheme_info, min_ref_frac, min_alt_frac,
                                  min_alt_freq)
    logger.info("Populating rule sets")
    valid_uids = list(assoc_data['kmer_rules'].keys())
    if only_update:
        rules = scheme_info['genotype_rule_sets']
        kmer_rules = {}
        for uid in assoc_data['kmer_rules']:
            kmer_rules[uid] = {'positive_genotypes':[],'partial_genotypes':[]}
        for genotype in rules:
            if genotype not in genotype_counts:
                for uid in rules[genotype]['positive_uids']:
                    kmer_rules[uid]['positive_genotypes'].append(genotype)
                for uid in rules[genotype]['partial_uids']:
                    kmer_rules[uid]['partial_genotypes'].append(genotype)
            else:
                for uid in assoc_data['kmer_rules']:
                    if genotype in assoc_data['kmer_rules'][uid]['positive_uids']:
                        kmer_rules[uid]['positive_genotypes'].append(genotype)
                    elif genotype in assoc_data['kmer_rules'][uid]['partial_genotypes']:
                        kmer_rules[uid]['partial_genotypes'].append(genotype)
        assoc_data['kmer_rules'] = kmer_rules

    else:
        logger.info("Filtering scheme sites which are present in less than {} samples".format(max_missing_count))
        valid_uids = list( set(valid_uids) - set(filter_missing_sites(kdata['kmer_counts'], scheme_info, max_missing_count)))
        if num_genotypes > 1:
            assoc_data['kmer_rules'] = blank_invalid_rules(assoc_data['kmer_rules'], num_genotypes, scheme_info)
        invalid_uids = set(list(scheme_info['uid_to_state'].keys())) - set(valid_uids)
        for uid in invalid_uids:
            assoc_data['kmer_rules'][uid] = {
            'positive_genotypes': [],
            'partial_genotypes':[]
        }

    if no_mixed:
        valid_uids = list(set(valid_uids) - set(list(kdata['mixed_sites'].keys())))

    invalid_uids = list(set(list(scheme_info['uid_to_state'].keys())) - set(valid_uids))

    out_str = []
    for uid in invalid_uids:
        row = []
        row.append(str(uid))
        row.append(str(scheme_info['uid_to_mutation'][uid]))
        row.append(str(scheme_info['uid_to_dna_name'][uid]))
        row.append(str(scheme_info['uid_to_aa_name'][uid]))
        row.append(str(scheme_info['uid_to_gene_feature'][uid]))
        row.append(str(scheme_info['uid_to_state'][uid]))
        row.append(str(scheme_info['uid_to_kseq'][uid]))
        out_str.append("{}".format("\t".join(row)))
    out_fh = open(os.path.join(outdir, "{}-removed-kmers.txt".format(prefix)), 'w')
    out_fh.write("{}".format("\n".join(out_str)))
    out_fh.close()
    stime = time.time()
    update_scheme(scheme_file, scheme_outfile, valid_uids, assoc_data['kmer_rules'], assoc_data['entropy'])
    print(time.time() - stime)

