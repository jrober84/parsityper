#!/usr/bin/python

from argparse import (ArgumentParser)
import logging, os, re, datetime, time
import pandas as pd
from parsityper.helpers import init_console_logger,read_fasta, read_samples
import multiprocessing as mp
from multiprocessing import Pool
from parsityper.scheme import parseScheme, constructSchemeLookups, SCHEME_HEADER
from parsityper.kmerSearch.kmerSearch import init_automaton_dict,perform_kmerSearch_fasta,perform_kmerSearch_fastq,process_kmer_results
from parsityper.version import __version__

if mp.get_start_method(allow_none=True) != 'spawn':
        mp.set_start_method('spawn', force=True)

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Parsityper genotype tuning')
    parser.add_argument('--input', type=str, required=True,
                        help='TSV file of sample_id,genotype,file_1,file_2')
    parser.add_argument('--profile', type=str, required=False,
                        help='Pre-computed kmer result profile for samples')
    parser.add_argument('--scheme', type=str, required=True,
                        help='TSV formated kmer scheme')
    parser.add_argument('--outdir', type=str, required=True,
                        help='output directory')
    parser.add_argument('--prefix', type=str, required=False,
                        help='output file prefix',default='parsityper')
    parser.add_argument('--update', type=str, required=False,
                        help='Do not blank existing rules, only update rules based on supplied samples')
    parser.add_argument('--min_members', type=int, required=False,
                        help='minimum number of representitives per genotype to be included',default=0)
    parser.add_argument('--min_cov', type=int, required=False,
                        help='Absolute minimum kmer coverage for fastq read detection default=auto determine coverage',default=50)
    parser.add_argument('--min_cov_frac', type=float, required=False,
                        help='Minimum percentage of total pool required for k-mer detection range 0 - 1.0 (default=0.05)', default=0.05)
    parser.add_argument('--min_partial_frac', type=float, required=False,
                        help='Minimum fraction of isolates positive for mutation for it to be partial 0 - 1.0 (default=0.1)', default=0.1)
    parser.add_argument('--min_alt_frac', type=float, required=False,
                        help='Minimum fraction of isolates positive for mutation for it to be positive 0 - 1.0 (default=0.1)', default=0.95)
    parser.add_argument('--min_ref_frac', type=float, required=False,
                        help='Minimum fraction of isolates positive for mutation for it to be positive 0 - 1.0 (default=0.1)', default=0.95)
    parser.add_argument('--min_positive_freq', type=int, required=False,
                        help='Minimum number of isolates positive for mutation for it to be valid for the scheme (default=1)', default=0)
    parser.add_argument('--max_frac_missing', type=float, required=False,
                        help='Maximum number of sequences allowed to be missing a mutation for it to stay included', default=0.25)
    parser.add_argument('--max_conflict', type=float, required=False,
                        help='Maximum percentage of genotype samples before droping as a requirement range 0 - 1.0 (default=0.05)', default=0.05)
    parser.add_argument('--n_threads', type=int, required=False,
                        help='output directory',default=1)
    parser.add_argument('-V', '--version', action='version', version='%(prog)s {}'.format(__version__))

    return parser.parse_args()

def summarize_samples(kmer_results,scheme_info,min_freq):
    scheme_global_tracker = {
        'mutation_freqs': {},
        'mutation_alt_count':{},
        'mutation_ave_cov': {},
        'mutation_sites_missing':{},
        'kmer_freqs': {},
        'kmer_missing': {},
    }
    num_samples = len(kmer_results)
    # Identify presence and absence of scheme mutations
    for sampleID in kmer_results:
        mutation_freq = kmer_results[sampleID]['total_mutation_key_freq']
        for mutation_key in mutation_freq:
            if not mutation_key in scheme_global_tracker['mutation_freqs']:
                scheme_global_tracker['mutation_freqs'][mutation_key] = 0
                scheme_global_tracker['mutation_sites_missing'][mutation_key] = []
                scheme_global_tracker['mutation_alt_count'][mutation_key] = 0
                scheme_global_tracker['mutation_ave_cov'][mutation_key] = 0

            scheme_global_tracker['mutation_freqs'][mutation_key] += mutation_freq[mutation_key]
            if mutation_freq[mutation_key] < min_freq:
                scheme_global_tracker['mutation_sites_missing'][mutation_key].append(sampleID)

        kmer_freq = kmer_results[sampleID]['raw_kmer_freq']
        for uid in kmer_freq:
            if not uid in scheme_global_tracker['kmer_freqs']:
                scheme_global_tracker['kmer_freqs'][uid] = 0

            freq = kmer_freq[uid]
            state = scheme_info['uid_to_state'][uid]
            mutation_key = scheme_info['uid_to_mutation'][uid]
            if freq < min_freq:
                if uid not in scheme_global_tracker['kmer_missing']:
                    scheme_global_tracker['kmer_missing'][uid] = []
                scheme_global_tracker['kmer_missing'][uid].append(sampleID)
            else:
                scheme_global_tracker['kmer_freqs'][uid] += 1
                if state == 'alt':
                    scheme_global_tracker['mutation_alt_count'][mutation_key] += 1

    for mutation_key in scheme_global_tracker['mutation_ave_cov']:
        ave_cov = 0
        total_freq = scheme_global_tracker['mutation_freqs'][mutation_key]
        total_samples = num_samples - len(scheme_global_tracker['mutation_sites_missing'][mutation_key])
        if total_samples > 0:
            ave_cov = total_freq / total_samples
        scheme_global_tracker['mutation_ave_cov'] = ave_cov

    return scheme_global_tracker

def identify_missing_kmers(genotype_rules,detected_scheme_kmers):
    return list(set(genotype_rules['positive_uids']) - set(detected_scheme_kmers))

def identify_conflict_kmers(genotype_rules,detected_scheme_kmers):
    return list(set(genotype_rules['positive_uids']) - set(detected_scheme_kmers))

def summarizeConflicts(sampleManifest,kmer_results,scheme_info,nthreads=1):
    genotype_conflicts = {}
    for sampleID in sampleManifest:
        genotype = sampleManifest[sampleID]['reported_genotype']
        if not sampleID in kmer_results:
            continue
        if not genotype in genotype_conflicts:
            genotype_conflicts[genotype] = {}
        if genotype in scheme_info['genotype_rule_sets']:
            genotype_rules = scheme_info['genotype_rule_sets'][genotype]
        else:
            genotype_rules = {'positive_uids':[],'positive_ref':[],'positive_alt':[]}
            scheme_info['genotype_rule_sets'][genotype] = genotype_rules
        detected_scheme_kmers = kmer_results[sampleID]['detected_scheme_kmers']
        sampleManifest[sampleID]['detected_scheme_kmers'] = detected_scheme_kmers
        sampleManifest[sampleID]['num_detected_scheme_kmers'] = len(detected_scheme_kmers)
        conflict_uids = identify_conflict_kmers(genotype_rules, detected_scheme_kmers)
        sampleManifest[sampleID]['conflict_kmers'] = conflict_uids
        for uid in conflict_uids:
            if not uid in genotype_conflicts[genotype]:
                genotype_conflicts[genotype][uid] = 0
            genotype_conflicts[genotype][uid]+=1
    return genotype_conflicts

def process_rawSamples(sampleManifest,scheme_info,nthreads,min_cov,report_sample_kmer_profiles):
    #Identify kmers in each sample
    if nthreads > 1:
        pool = Pool(processes=nthreads)
    logging.info("Performing kmer searching on {} samples with {} threads".format(len(sampleManifest) ,nthreads))

    # Init Ahocorasak automation objects
    logging.info("Initializing aho-corasick automation")
    aho = {'scheme': init_automaton_dict(scheme_info['uid_to_kseq'])}
    kmer_results = {}
    for sampleID in sampleManifest:
        seq_files = sampleManifest[sampleID]['raw_seq_files']
        fileType = sampleManifest[sampleID]['file_type']

        if nthreads == 1:
            if fileType == 'fastq':
                kmer_results[sampleID] = perform_kmerSearch_fastq(scheme_info['uid_to_kseq'],
                                                                  scheme_info['kseq_to_uids'],
                                                                  aho['scheme'], seq_files)
            else:
                kmer_results[sampleID] = perform_kmerSearch_fasta(scheme_info['uid_to_kseq'],
                                                                  scheme_info['kseq_to_uids'],
                                                                  aho['scheme'], read_fasta(seq_files[0]), min_cov)
        else:
            if fileType == 'fastq':
                kmer_results[sampleID] = pool.apply_async(perform_kmerSearch_fastq, (
                    scheme_info['uid_to_kseq'], scheme_info['kseq_to_uids'], aho['scheme'], seq_files))
            else:
                kmer_results[sampleID] = pool.apply_async(perform_kmerSearch_fasta,
                                                          (scheme_info['uid_to_kseq'], scheme_info['kseq_to_uids'],
                                                           aho['scheme'], read_fasta(seq_files[0]), min_cov))

    # Extract results in multithreaded mode
    if nthreads > 1:
        pool.close()
        pool.join()
        for sampleID in kmer_results:
            kmer_results[sampleID] = kmer_results[sampleID].get()

    logging.info("Kmer searching complete")
    logging.info("Creating kmer profiles")
    kmer_results_df = pd.DataFrame.from_dict(kmer_results,orient='index').transpose()
    logging.info("Writting kmer profiles")
    kmer_results_df.to_csv(report_sample_kmer_profiles, sep="\t", header=True)
    return kmer_results

def init_sampleManifest(samples,scheme_name,scheme_info,analysis_date):
    sampleManifest = {}
    for fileType in samples:
        for i in range(0,len(samples[fileType])):
            row = samples[fileType][i]
            sampleID = row['sample_id']
            seq_files = [row['file_1']]
            if row['file_2'] != '':
                seq_files.append(row['file_2'])

            sampleManifest[sampleID] = {
                'sample_id': sampleID,
                'scheme': scheme_name,
                'analysis_date': analysis_date,
                'reported_genotype': row['genotype'],
                'file_type': fileType,
                'num_seq_files': len(seq_files),
                'raw_seq_files': seq_files,
                'compatible_genotypes': [],
                'total_scheme_kmers': len(scheme_info['uid_to_kseq']),
                'detected_scheme_kmers': [],
                'num_detected_scheme_kmers': 0,
                'conflict_kmers': [],
                'ave_scheme_kmers_freq': 0,
                'total_scheme_mutations': len(scheme_info['mutation_to_uid']),
                'detected_scheme_mutations': 0,
                'detected_scheme_mixed_mutations': 0,
                'primary_genotype': '',
                'primary_genotype_frac': '',
            }
    return sampleManifest

def updateSchemeInfo(scheme_info,valid_mutations,valid_uids,valid_genotypes):
    scheme_info['genotypes'] = valid_genotypes
    profiles = {
        'max_variant_positions': 0,
        'gene_features': [],
        'genotypes': [],
        'genotype_rule_sets': scheme_info['genotype_rule_sets'],
        'kmer_to_genotypes': {},
        'uid_to_state': {},
        'uid_to_kseq': {},
        'kseq_to_uids': {},
        'uid_to_mutation': {},
        'uid_to_dna_name': {},
        'uid_to_aa_name': {},
        'uid_to_gene_feature': {},
        'mutation_to_uid': {},
        'kmer_profiles': {},
        'mutation_profiles': {},
        'min_kmer_len': 100000,
        'max_kmer_len': 0
    }
    min_kmer_len = 10000
    max_kmer_len = 0
    new_uid = 0
    for uid in valid_uids:
        if not uid in scheme_info['uid_to_kseq']:
            continue
        if not uid in scheme_info['uid_to_mutation']:
            continue
        mutation = scheme_info['uid_to_mutation'][uid]
        if not mutation in valid_mutations:
            continue

        kseq = scheme_info['uid_to_kseq'][uid]
        state = scheme_info['uid_to_state'][uid]
        dna_name= scheme_info['uid_to_dna_name'][uid]
        aa_name = scheme_info['uid_to_aa_name'][uid]
        feature = scheme_info['uid_to_gene_feature'][uid]

        if not kseq in profiles['kseq_to_uids']:
            profiles['kseq_to_uids'][kseq] = []

        if not mutation in profiles['mutation_to_uid']:
            profiles['mutation_to_uid'][mutation] = []

        profiles['mutation_to_uid'][mutation].append(new_uid)
        profiles['kseq_to_uids'][kseq].append(new_uid)
        profiles['uid_to_state'][new_uid] = state
        profiles['uid_to_kseq'][new_uid] = kseq
        profiles['uid_to_dna_name'][new_uid] = dna_name

        profiles['uid_to_dna_name'][new_uid] = dna_name
        profiles['uid_to_aa_name'][new_uid] = aa_name
        profiles['uid_to_gene_feature'][new_uid] = feature

        klen= len(kseq)
        if klen < min_kmer_len:
            min_kmer_len = klen
        if klen > max_kmer_len:
            max_kmer_len = klen
        new_uid+=1
    profiles['min_kmer_len'] = min_kmer_len
    profiles['max_kmer_len'] = max_kmer_len

    for genotype in profiles['genotype_rule_sets']:
        profiles['genotype_rule_sets'][genotype]['positive_uids'] = list(
            set(profiles['genotype_rule_sets'][genotype]['positive_uids']) & set(valid_uids))
        profiles['genotype_rule_sets'][genotype]['positive_ref'] = list(
            set(profiles['genotype_rule_sets'][genotype]['positive_ref']) & set(valid_uids))
        profiles['genotype_rule_sets'][genotype]['positive_alt'] = list(
            set(profiles['genotype_rule_sets'][genotype]['positive_alt']) & set(valid_uids))



    mutations = list(profiles['mutation_to_uid'].keys())
    num_mutations = len(mutations)

    # populate the profiles
    for i in range(0, len(valid_genotypes)):
        genotype = valid_genotypes[i]
        profiles['mutation_profiles'][genotype] = [0.5] * num_mutations

        for k in range(0, num_mutations):
            uids = profiles['mutation_to_uid']
            ref_count = len(set(uids)  & set(profiles['genotype_rule_sets'][genotype]['positive_ref']))
            alt_count = len(set(uids) & set(profiles['genotype_rule_sets'][genotype]['positive_alt']))
            value = -1
            if ref_count > 0 and alt_count > 0:
                value = 0.5
            elif ref_count > 0:
                value = 0
            elif alt_count > 0:
                value = 1
            profiles['mutation_profiles'][genotype][k] = value

    return profiles

def updateScheme(scheme_file,scheme_info,outfile):
    df = pd.read_csv(scheme_file, sep="\t", header=0)
    fh = open(outfile,'w')
    fh.write("{}\n".format("\t".join(SCHEME_HEADER)))
    columns = df.columns.tolist()
    num_fields = len(SCHEME_HEADER)
    new_uid_key = 0
    rules = {}

    for genotype in scheme_info['genotype_rule_sets']:
        uids = scheme_info['genotype_rule_sets'][genotype]['positive_uids']
        for uid in uids:
            if not uid in rules:
                rules[uid] = []
            state = scheme_info['uid_to_state'][uid]
            if state == 'alt':
                rules[uid].append(genotype)

    for index,row in df.iterrows():
        uid = row['key']
        mutation_key = row['mutation_key']
        if not mutation_key in scheme_info['mutation_to_uid']:
            continue
        seq = row['unalign_kseq']
        if not seq in scheme_info['kseq_to_uids']:
            continue
        entry = {}
        for field in SCHEME_HEADER:
            if field in columns:
                value = str(row[field])
                if value == 'nan':
                    value = ''
                entry[field] = value
            else:
                entry[field] = ''
        entry['key'] = new_uid_key
        positive_genotypes = []
        if uid in rules:
            positive_genotypes = rules[uid]
            #print("{}\t{}".format(uid,rules[uid]))


        entry['positive_genotypes'] = ','.join([str(x) for x in positive_genotypes])
        entry['seq_ids'] = ''
        entry['partial_genotypes'] = ''

        record = []
        for i in range(0,num_fields):
            record.append(entry[SCHEME_HEADER[i]])
        fh.write("{}\n".format("\t".join([str(x) for x in record])))
        new_uid_key+=1
    fh.close()

def compare_sample_to_genotype(data,genotype,genotype_results,scheme_info):
    uids = set(list(scheme_info['uid_to_mutation'].keys()))
    valid_uids = set(data['valid_uids'])
    mixed_sites = list(data['mixed_sites'])
    missing_sites = list(uids - valid_uids )
    exclude_sites = set(missing_sites + mixed_sites)
    mixed_sites = set(mixed_sites)
    missing_sites = set(missing_sites)
    detected_scheme_kmers = set(data['detected_scheme_kmers'])
    geno_rules = scheme_info['genotype_rule_sets']

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

    results = {
        'dist':0,
        'missing_sites':missing_sites,
        'matched_kmers': [],
        'matched_alt_kmers': [],
        'mismatch_kmers':[],
        'mismatch_alt_kmers':[],
        'is_compatible':True
    }

    dist = 0
    uids = valid_uids & set(geno_rules[genotype]['positive_uids']) - exclude_sites
    mismatches = uids - detected_scheme_kmers
    results['mismatch_kmers'] = mismatches
    matched = list(detected_scheme_kmers & (uids - mismatches))
    results['matched_kmers'] = matched
    positive_alt_match = list(set(geno_rules[genotype]['positive_alt']) & set(matched ))
    results['matched_alt_kmers'] = positive_alt_match
    positive_alt_mismatch = list(set(geno_rules[genotype]['positive_alt']) & set(mismatches))
    results['mismatch_alt_kmers'] = positive_alt_mismatch

    genotype_results[genotype]['matched_pos_kmers'] = matched
    if len(uids) > 0:
        dist =  len(mismatches) / len(uids)
    results['dist'] = dist

    if len(mismatches) > 0:
        results['is_compatible'] = False

    return results


def run():
    cmd_args = parse_args()
    logger = init_console_logger(2)

    #input parameters
    prefix = cmd_args.prefix
    min_cov = int(cmd_args.min_cov)
    min_cov_frac = float(cmd_args.min_cov_frac)
    min_members = int(cmd_args.min_members)
    only_update = cmd_args.update
    input = cmd_args.input
    profile = cmd_args.profile
    min_alt_frac = cmd_args.min_alt_frac
    min_ref_frac = cmd_args.min_ref_frac
    min_partial_frac = cmd_args.min_partial_frac
    min_pos_freq = cmd_args.min_positive_freq
    max_conflict = cmd_args.max_conflict

    scheme_file = cmd_args.scheme
    max_frac_missing = cmd_args.max_frac_missing
    nthreads = cmd_args.n_threads
    outdir = cmd_args.outdir
    logger = init_console_logger(2)

    #result files
    report_sample_kmer_profiles = os.path.join(outdir,"sample_kmer.profiles.txt")
    scheme_outfile = os.path.join(outdir,"{}-scheme.txt".format(prefix))

    # initialize analysis directory
    if not os.path.isdir(outdir):
        logger.info("Creating analysis results directory {}".format(outdir))
        os.mkdir(outdir, 0o755)
    else:
        logger.info("Results directory {} already exits, will overwrite any results files here".format(outdir))

    scheme = parseScheme(scheme_file)

    logger.info("Initializing scheme data structure from {}".format(scheme_file))
    scheme_info = constructSchemeLookups(scheme)

    scheme_name = os.path.basename(scheme_file)
    analysis_date = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    logger.info("Initializing sample manifest from {}".format(input))
    samples = read_samples(input)

    sampleManifest = init_sampleManifest(samples,scheme_name,scheme_info,analysis_date)

    num_samples = len(sampleManifest )
    if profile == None:
        kmer_results = process_rawSamples(sampleManifest,scheme_info,nthreads,min_cov,report_sample_kmer_profiles)
    else:
        kmer_results_df = pd.read_csv(profile,sep="\t",header=True)
        sample_list = kmer_results_df.columns.tolist()
        #TODO implement profile check

    #process kmer results
    logger.info("Processing kmer results")
    kmer_results = process_kmer_results(scheme_info, kmer_results, min_cov, min_cov_frac)
    logger.info("Summarizing kmer results by genotype")
    scheme_kmer_result_summary = {'dataset':summarize_samples(kmer_results,scheme_info,min_cov),'genotypes':{}}


    invalid_mutation_keys = []

    for mutation_key in scheme_kmer_result_summary['dataset']['mutation_sites_missing']:
        frac_missing = len(scheme_kmer_result_summary['dataset']['mutation_sites_missing'][mutation_key]) / num_samples
        if frac_missing > max_frac_missing:
            invalid_mutation_keys.append(mutation_key)
    logger.info("Found {} mutations where it is present in less than {}% samples".format(len(invalid_mutation_keys),max_frac_missing))

    #get counts of genotypes
    genotypeCounts = {}
    genotypeMap = {}
    for sampleID in sampleManifest:
        genotype = sampleManifest[sampleID]['reported_genotype']
        if not genotype in genotypeCounts:
            genotypeCounts[genotype] = 0
            genotypeMap[genotype] = []
        genotypeCounts[genotype] += 1
        genotypeMap[genotype].append(sampleID)

    #get summary of each genotype
    for genotype in genotypeMap:
        subset = {}
        for sample_id in genotypeMap[genotype]:
            subset[sample_id] = kmer_results[sample_id]
        scheme_kmer_result_summary['genotypes'][genotype] = summarize_samples(subset,scheme_info,min_cov)

    genotype_conflicts_pre = summarizeConflicts(sampleManifest, kmer_results, scheme_info, nthreads)
    new_rules = {}
    for genotype in genotypeCounts:
        kmer_counts = scheme_kmer_result_summary['genotypes'][genotype]['kmer_freqs']
        new_rules[genotype] = {'positive_uids':[],'positive_ref':[],'positive_alt':[]}
        gCount  = genotypeCounts[genotype]
        for uid in kmer_counts:
            state = scheme_info['uid_to_state'][uid]

            freq = kmer_counts[uid]

            perc_present = 0
            if gCount > 0:
                perc_present = freq /  gCount

            if perc_present >= min_alt_frac :

                if state == 'ref':
                    new_rules[genotype]['positive_ref'].append(uid)
                else:
                    new_rules[genotype]['positive_uids'].append(uid)
                    new_rules[genotype]['positive_alt'].append(uid)
                    print("{}\t{}\t{}".format(uid,state,perc_present))
        scheme_info['genotype_rule_sets'][genotype] = new_rules[genotype]

    print(scheme_info['genotype_rule_sets']['A.1']['positive_alt'])
    genotype_conflicts_post = summarizeConflicts(sampleManifest, kmer_results, scheme_info, nthreads)


    #get valid genotypes
    valid_genotypes = []
    invalid_mutation_keys = []
    if not only_update:
        for genotype in genotypeCounts:
            count = genotypeCounts[genotype]
            if count < min_members:
                continue
            valid_genotypes.append(genotype)
        valid_genotypes.sort()

        invalid_genotypes = list(set(scheme_info['genotypes']) - set(valid_genotypes))
        logger.info("Found {} genotypes with too few members, these will be removed from the scheme".format(len(invalid_genotypes)))

        #filter out invalid genotypes from scheme
        scheme_info['genotypes'] = valid_genotypes
        for genotype in invalid_genotypes:
            if genotype in scheme_info['genotype_rule_sets']:
                del(scheme_info['genotype_rule_sets'][genotype])


        logger.info("Removing kmers which are present in too few samples")
        #Get invalid kmers
        invalid_kmers = []
        for uid in scheme_kmer_result_summary['dataset']['kmer_missing']:
            count_missing = len(scheme_kmer_result_summary['dataset']['kmer_missing'][uid])
            count_present = num_samples - count_missing
            if count_present < min_members:
                invalid_kmers.append(uid)

        logger.info("Flagged {} kmers for removal due to presence in too few samples".format(len(invalid_kmers)))
        logger.info("Identifying mutations which do not have both alt/ref states")

        # flag mutations which are now only one state due to kmer filtering
        for mutation_key in scheme_info['mutation_to_uid']:
            uids = scheme_info['mutation_to_uid'][mutation_key]
            is_ref_present = False
            is_alt_present = False
            for uid in uids:
                state = scheme_info['uid_to_state'][uid]
                if state == 'ref':
                    is_ref_present = True
                else:
                    is_alt_present = True
                if is_ref_present and is_alt_present:
                    break
            if not is_ref_present or not is_alt_present:
                invalid_mutation_keys.append(mutation_key)

        logger.info("Flagged {} mutations for removal where both alt/ref states are not present".format(len(invalid_mutation_keys)))
        invalid_mutation_keys = list(set(invalid_mutation_keys))

        #filter out invalid mutations from scheme
        for mutation_key in invalid_mutation_keys:
            if not mutation_key in scheme_info['mutation_to_uid']:
                continue
            uids = scheme_info['mutation_to_uid'][mutation_key]
            invalid_kmers.extend(uids)
            for uid in uids:
                fields = list(scheme_info.keys())
                for field in fields:
                    if 'uid' in field:
                        if uid in scheme_info[field]:
                            del(scheme_info[field][uid])
            del(scheme_info['mutation_to_uid'][mutation_key])

        #update the scheme info
        valid_mutations = list(set(list(scheme_info['mutation_to_uid'].keys())) - set(invalid_mutation_keys))
        logger.info("found {} valid mutations after filtering".format(len(valid_mutations)))
        valid_uids = []
        for mutation_key in valid_mutations:
            valid_uids.extend(scheme_info['mutation_to_uid'][mutation_key])
        valid_uids = list(set(valid_uids) - set(invalid_kmers))
        logger.info("found {} valid kmers after filtering".format(len(valid_uids)))

    logger.info("Writting updated scheme to {}".format(scheme_outfile))
    updateScheme(scheme_file, scheme_info, scheme_outfile)

