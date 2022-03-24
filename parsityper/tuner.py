#!/usr/bin/python
import logging
import multiprocessing as mp
import os
import sys
from argparse import (ArgumentParser)
from scipy.stats import entropy
from parsityper.helpers import init_console_logger, read_tsv
from parsityper.scheme import parseScheme, constructSchemeLookups
from parsityper.version import __version__

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Parsityper genotype tuning')
    parser.add_argument('--input_meta', type=str, required=True,
                        help='TSV file of sample_id,genotype')
    parser.add_argument('--input_profile', type=str, required=False,
                        help='Pre-computed kmer result profile for samples')
    parser.add_argument('--scheme', type=str, required=True,
                        help='TSV formated kmer scheme')
    parser.add_argument('--outdir', type=str, required=True,
                        help='output directory')
    parser.add_argument('--prefix', type=str, required=False,
                        help='output file prefix',default='parsityper')
    parser.add_argument('--update', required=False,
                        help='Do not blank existing rules, only update rules based on supplied samples',action='store_true')
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

def get_genotype_kmer_counts(sample_mapping,kmer_file,min_cov=1):
    genotypes = sorted(list(set(list(sample_mapping.values()))))
    genotype_kmer_counts = {}
    for genotype in genotypes:
        genotype_kmer_counts[genotype] = []
    uid = 0
    kmer_counts = {}
    kmer_FH = open(kmer_file,'r')
    samples = next(kmer_FH).rstrip().split("\t")[1:]
    sample_range = range(1,len(samples))
    for line in kmer_FH:
        if not uid in kmer_counts:
            kmer_counts[uid] = 0
        for genotype in genotypes:
            genotype_kmer_counts[genotype].append(0)
        row = line.rstrip().split(0)
        for i in sample_range:
            sample_id = samples[i]
            genotype = sample_mapping[sample_id]
            value = row[i]
            if value > min_cov:
                genotype_kmer_counts[genotype][uid]+=1
                kmer_counts[uid]+=1

        uid+=1
    return {'genotype_kmer_counts':genotype_kmer_counts,'kmer_counts':kmer_counts}

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
    for genotype in genotype_kmer_counts:
        genotype_count = genotype_counts[genotype]
        rules[genotype] = {'positive_uids':[],
                           'positive_ref':[],
                           'positive_alt':[],
                           'partial_uids':[],
                           'partial_ref':[],
                           'partial_alt':[]}
        for uid,value in enumerate(genotype_kmer_counts[genotype]):
            if value == 0:
                continue
            state = scheme_info['uid_to_state'][uid]
            frac = value / genotype_count
            if frac < min_par_frac:
                continue
            if state == 'ref' and frac >= min_ref_frac:
                rules[genotype]['positive_uids'].append(uid)
                rules[genotype]['positive_ref'].append(uid)
            elif state == 'alt' and frac >= min_alt_freq:
                rules[genotype]['positive_uids'].append(uid)
                rules[genotype]['positive_alt'].append(uid)
            elif state == 'ref' and frac >= min_par_frac:
                rules[genotype]['partial_uids'].append(uid)
                rules[genotype]['positive_ref'].append(uid)
            elif state == 'alt' and frac >= min_par_frac:
                rules[genotype]['partial_uids'].append(uid)
                rules[genotype]['partial_alt'].append(uid)

        sEntropy[uid] = calc_shanon_entropy(list(genotype_kmer_counts[genotype].values()))

    kmer_rules = {}
    for uid,value in enumerate(genotype_kmer_counts[genotype]):
        kmer_rules[uid] = {
            'positive_genotypes': [],
            'partial_genotypes':[]
        }

    for genotype in rules:
        for uid in rules[genotype]['positive_uids']:
            kmer_rules[uid]['positive_genotypes'].append(genotype)
        for uid in rules[genotype]['partial_uids']:
            kmer_rules[uid]['partial_genotypes'].append(genotype)

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

def filter_missing_sites(kmer_counts,scheme_info,max_missing_count):
    invalid_kmers = []
    for mutation_key in scheme_info['mutation_to_uid']:
        uids = scheme_info['mutation_to_uid'][mutation_key]
        count_present = 0
        for uid in uids:
             count_present+= kmer_counts[uid]
        if count_present > max_missing_count:
            invalid_kmers += uids

    return invalid_kmers


def update_scheme(input_scheme,output_scheme,valid_uids,kmer_geno_rules,kmer_entropies):
    in_FH = open(input_scheme,'r')
    out_FH = open(output_scheme, 'w')
    header = next(in_FH).rtrim().split('\t')
    out_FH.write("{}\n".format("\t".join(header)))
    uid_col = 0
    entropy_col = header.index('entropy')
    pos_col = header.index('positive_genotypes')
    par_col = header.index('partial_genotypes')
    uid = 0
    for line in in_FH:
        row = line.rtrim().split('\t')
        row_uid = int(row[uid_col])
        if not row_uid in valid_uids:
            row[uid_col] = uid
        row[entropy_col] = kmer_entropies[row_uid]
        row[pos_col] = ','.join(kmer_geno_rules[row_uid]['positive_genotypes'])
        row[par_col] = ','.join(kmer_geno_rules[row_uid]['partial_genotypes'])
        out_FH.write("{}\n".format("\t".join([str(x) for x in row])))
        uid+=1
    out_FH.close()
    in_FH.close()


def run():
    # input parameters
    cmd_args = parse_args()
    input_meta = cmd_args.input_meta
    input_profile = cmd_args.in_profile
    scheme_file = cmd_args.scheme
    prefix = cmd_args.prefix
    min_cov = int(cmd_args.min_cov)
    only_update = cmd_args.update
    min_alt_frac = cmd_args.min_alt_frac
    min_ref_frac = cmd_args.min_ref_frac
    min_alt_freq = cmd_args.min_positive_freq

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
    max_missing_count = int(max_frac_missing * len(sample_mapping))

    #Read scheme
    logger.info("Initializing scheme data structure from {}".format(scheme_file))
    if not os.path.isfile(scheme_file):
        logger.error("Error scheme file {} does not exist".format(scheme_file))
        sys.exit()
    scheme = parseScheme(scheme_file)
    scheme_info = constructSchemeLookups(scheme)

    #Read Kmer header
    if not os.path.isfile(input_profile):
        logger.error("Error profile {} does not exist".format(input_profile))
        sys.exit()

    kmer_FH = open(input_profile,'r')
    header = next(kmer_FH).rstrip().split("\t")
    kmer_FH.close()
    profile_samples = header[1:]

    #Confirm all samples in profile exist in metadata
    ovl = set(profile_samples) & set(list(sample_mapping.keys()))
    if len(ovl) != len(set(profile_samples)):
        logger.error("Sample id's in profile are not present in metadata {}".format(set(profile_samples) - set(list(sample_mapping.keys()))))
        sys.exit()

    logger.info("Found {} samples in profile {}".format(len(profile_samples),input_profile))

    kdata = get_genotype_kmer_counts(sample_mapping, input_profile, min_cov)

    genotype_counts = get_genotype_counts(sample_mapping)
    num_genotypes = len(genotype_counts)
    assoc_data = determine_genotype_kmer_assoc(kdata['genotype_kmer_counts'], genotype_counts, scheme_info, min_ref_frac, min_alt_frac,
                                  min_alt_freq)

    assoc_data['kmer_rules'] = blank_invalid_rules(assoc_data['kmer_rules'], num_genotypes, scheme_info)
    valid_uids = list(assoc_data['kmer_rules'].keys())
    if only_update:
        rules = scheme_info['genotype_rule_sets']
        kmer_rules = {}
        for uid in assoc_data['kmer_rules']:
            kmer_rules[uid] = {'positive_uids':[],'partial_genotypes':[]}
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
        valid_uids = list( set(valid_uids) - set(filter_missing_sites(kdata['kmer_counts'], scheme_info, max_missing_count)))


    update_scheme(scheme_file, scheme_outfile, valid_uids, assoc_data['kmer_rules'], assoc_data['entropy'])

