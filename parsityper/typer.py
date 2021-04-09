#!/usr/bin/python
from argparse import (ArgumentParser, FileType)
import logging, os, sys, re, collections, operator, math, shutil, datetime
from collections import Counter
import pandas as pd
from parsityper.helpers import validate_args, init_console_logger, read_tsv, scheme_to_biohansel_fasta,filter_biohansel_kmer_df, \
    init_kmer_targets,generate_biohansel_kmer_names, get_scheme_template

from parsityper.bio_hansel import bio_hansel
import statistics
from parsityper.constants import PRIMER_SCHEMES, TYPING_SCHEMES

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='SARS-COV-2 Kmer-analysis')
    parser.add_argument('--data_dir', type=str, required=False,
                        help='directory of fasta/fastq files')
    parser.add_argument('--mode', type=str, required=True,
                        help='Operate in batch or single mode. Incompatible with input options other than data_dir')
    parser.add_argument('--type', type=str, required=True,
                        help='Treat data as single or multi source')
    parser.add_argument('--se', type=str, required=False,
                        help='single-end fastq read')
    parser.add_argument('--R1', type=str, required=False,
                        help='paired-end fwd fastq read')
    parser.add_argument('--R2', type=str, required=False,
                        help='paired-end rev fastq read')
    parser.add_argument('--scheme', type=str, required=False,
                        help='TSV formated kmer scheme',default='SARS-COV-2')
    parser.add_argument('--primers', type=str, required=False,
                        help='TSV formated primer file',default='arctic')
    parser.add_argument('--outdir', type=str, required=True,
                        help='output directory')
    parser.add_argument('--prefix', type=str, required=False,
                        help='output file prefix',default='SARS-COV-2_kmer_analysis')
    parser.add_argument('--min_cov', type=int, required=False,
                        help='Absolute minimum kmer coverage for fastq read detection default=auto determine coverage',default=50)
    parser.add_argument('--min_cov_frac', type=float, required=False,
                        help='Minimum percentage of total pool required for k-mer detection range 0 - 1.0 (default=0.05)', default=0.05)
    parser.add_argument('--max_mixed_sites', type=int, required=False,
                        help='Maximum number of sites allowed to have both kmer states', default=5)
    parser.add_argument('--max_missing_sites', type=int, required=False,
                        help='Maximum number of sites allowed to be missing', default=500)
    parser.add_argument('--sample_dist_cutoff', type=float, required=False,
                        help='Maximum profile distance for including sample in a comparison', default=0.05)
    parser.add_argument('--genotype_dist_cutoff', type=float, required=False,
                        help='Maximum profile distance for considering a genotype result valid', default=0.05)
    parser.add_argument('--n_threads', type=str, required=False,
                        help='output directory',default=1)
    parser.add_argument('--strain_profiles', type=str, required=False,
                        help='Kmer profile scheme to compare against')
    parser.add_argument('--genotype_profiles', type=str, required=False,
                        help='Kmer profile scheme to compare against')
    parser.add_argument('--no_template_control', type=str, required=False, help='Fasta/Fastq formatted no template control data')
    parser.add_argument('--positive_control', type=str, required=False,
                        help='Fasta/Fastq formatted positive control data')
    parser.add_argument('--list', type=str, required=False,
                        help='list in-built primer and typing schemes')
    return parser.parse_args()


def generate_sample_kmer_profile(biohansel_df,scheme_kmer_target_keys,min_cov=20):
    '''
    Accepts biohansel kmer pandas dataframe and sums the frequency of positive and negative kmers by sample if they meet the minimum threshold
    :param biohansel_df: pandas biohansel kmer df
    :param scheme_kmer_target_keys: list of kmer targets
    :param min_cov: integer of minimum frequency
    :return: dict of samples and counts of positive and negative kmers
    '''
    samples = {}
    for row in biohansel_df.itertuples():
        sample = row.sample
        is_pos_kmer = row.is_pos_kmer
        target = row.kmername.replace('negative','').split('-')[0]

        freq = row.freq
        if not sample in samples:
            samples[sample] = {}
            samples[sample]['counts'] = get_scheme_template(scheme_kmer_target_keys,{'positive':0,'negative':0})
        if not target in samples[sample]['counts']:
            continue
        if freq > min_cov:
            if is_pos_kmer:
                samples[sample]['counts'][target]['positive'] += freq
            else:
                samples[sample]['counts'][target]['negative'] += freq
    #print(samples[sample])


    return samples

def calc_kmer_ratio(biohansel_df,scheme_kmer_target_keys,min_cov=20):
    '''
    Accepts biohansel kmer pandas dataframe and sums the frequency of positive and negative kmers by sample if they meet the minimum threshold
    and then determines the ratio of the positive and negative kmers
    :param biohansel_df: pandas biohansel kmer df
    :param scheme_kmer_target_keys: list of kmer targets
    :param min_cov: integer of minimum frequency
    :return: dict of samples and counts of positive and negative kmers
    '''
    samples = generate_sample_kmer_profile(biohansel_df,scheme_kmer_target_keys,min_cov=20)
    for sample in samples:
        samples[sample]['ratios'] = {}
        for target in samples[sample]['counts']:
            total = sum(samples[sample]['counts'][target].values())
            if total == 0:
                ratio = -1
            else:
                ratio = samples[sample]['counts'][target]['positive'] / total
            samples[sample]['ratios'][target] = ratio
    return samples


def calc_mixed_sites(sample_data,min_cov_frac):
    '''
    Using sample data dictionary of kmer targets determines which sites are mixed by as determined by the cutoff value
    :param sample_data: dict
    :param min_cov_frac: float (0-1)
    :return: sample data dict with additional fields added
    '''
    for sample in sample_data:
        sample_data[sample]['statistics'] = {}
        sample_data[sample]['statistics']['num_kmers_targets_present'] = 0
        sample_data[sample]['statistics']['num_mixed_kmers_targets'] = 0
        sample_data[sample]['statistics']['mixed_kmers_targets'] = []
        total_freqs = []
        mixed_ratios = []
        for target in sample_data[sample]['ratios']:
            ratio = sample_data[sample]['ratios'][target]
            if ratio == -1:
                continue
            total_freqs.append(sum(sample_data[sample]['counts'][target].values()))
            sample_data[sample]['statistics']['num_kmers_targets_present']+=1
            if ratio == 1 or ratio == 0:
                continue

            if ratio >= min_cov_frac or ratio <= (1-min_cov_frac):
                sample_data[sample]['statistics']['num_mixed_kmers_targets']+=1
                sample_data[sample]['statistics']['mixed_kmers_targets'].append(target)
                mixed_ratios.append(ratio)
        if len(total_freqs) > 0:
            sample_data[sample]['statistics']['average_target_coverage'] = sum(total_freqs) / len(total_freqs)
            sample_data[sample]['statistics']['stdev_target_coverage'] = statistics.pstdev(total_freqs)
        else:
            sample_data[sample]['statistics']['average_target_coverage'] = 0
            sample_data[sample]['statistics']['stdev_target_coverage'] = 0

        if len(mixed_ratios) >0:
            sample_data[sample]['statistics']['average_mixed_target_ratio'] = sum(mixed_ratios) / len(mixed_ratios)
            sample_data[sample]['statistics']['stdev_mixed_target_ratio'] = statistics.pstdev(mixed_ratios)
        else:
            sample_data[sample]['statistics']['average_mixed_target_ratio'] = 0
            sample_data[sample]['statistics']['stdev_mixed_target_ratio'] = 0



    return sample_data

def get_mixed_kmer_targets(sample_data,min_cov_frac):
    '''
    Accepts sample kmer profile dict and determines for each target in the scheme the number of times that target was found to be mixed
    and which samples were involved
    :param sample_data: sample data dict
    :param min_cov_frac: float (0-1)
    :return:
    '''
    targets = {}
    for sample in sample_data:
        for target in sample_data[sample]['ratios']:
            if target not in targets:
                targets[target] = {
                    'count_samples_present':0,
                    'count_samples_absent': 0,
                    'count_mixed_samples': 0,
                    'mixed_samples_ids':[]
                }
            ratio = sample_data[sample]['ratios'][target]
            if ratio == -1:
                targets[target]['count_samples_absent']+=1
            else:
                targets[target]['count_samples_present'] += 1
                if (ratio != 0 and ratio !=1) and (ratio >= min_cov_frac or ratio <= (1-min_cov_frac)):
                    targets[target]['count_mixed_samples']+=1
                    targets[target]['mixed_samples_ids'].append(sample)
    return targets

def subtract_contamination_kmers(sample_data,kmers,min_cov_frac):
    '''

    :param sample_data:
    :param kmers:
    :param min_cov_frac:
    :return:
    '''

    for sample in sample_data:
        for kmer in kmers:
            if 'negative' in kmer:
                id = 'negative'
            else:
                id = 'positive'
            target = kmer.replace('negative','').split('-')[0]
            frequency = kmers[kmer]
            #The level of contamination in this sample is below what would be called at a sub-consensus level so no need to change
            #anything
            if sample_data[sample]['ratios'][target] < min_cov_frac or sample_data[sample]['ratios'][target] > 1 - min_cov_frac:
                continue
            sample_data[sample]['counts'][target][id] -= frequency
            if sample_data[sample]['counts'][target][id] < 0:
                sample_data[sample]['counts'][target][id] = 0
            total = sum(sample_data[sample]['counts'][target].values())
            if total == 0:
                sample_data[sample]['ratios'][target] = -1
            else:
                sample_data[sample]['ratios'][target] = sample_data[sample]['counts'][target]['positive']/total
    return sample_data

def identify_compatible_types(scheme_df,sample_data,min_cov_frac,detection_limit=100):
    '''
    Using the genotype information from the scheme dataframe, genotypes are determined to be compatible or incompatible
    with the available data provided that the site has sufficient kmer frequency
    :param scheme_df:
    :param sample_data:
    :param min_cov_frac:
    :return:
    '''
    genotype_kmer_mapping = get_scheme_genotypes(scheme_df)
    genotypes = list(genotype_kmer_mapping.keys())
    for row in scheme_df.itertuples():
        if 'nan' in str(row.positive_seqs):
            continue
        positive_seqs = str(row.positive_seqs).split(',')
        partial_pos_seqs = str(row.partial_positive_seqs).split(',')



        target = str(row.key)

        for sample in sample_data:
            if not 'genotypes' in sample_data[sample]:
                sample_data[sample]['genotypes'] = {
                    'include': [],
                    'informative':[],
                    'exclude':[],
                    'unique':[],
                    'candidates':[],
                    'candidate_data':{}

                }
            if target not in sample_data[sample]['counts']:
                continue

            ratio = sample_data[sample]['ratios'][target]
            if ratio == -1:
                continue
            if ratio > 0:
                if len(positive_seqs) < len(genotypes) * 0.75:
                    sample_data[sample]['genotypes']['include'].extend(positive_seqs)

                if len(partial_pos_seqs) < len(genotypes) * 0.75:
                    sample_data[sample]['genotypes']['include'].extend(partial_pos_seqs)
                if len(positive_seqs) < len(genotypes) * 0.01 and len(partial_pos_seqs) < len(genotypes) * 0.01:
                    sample_data[sample]['genotypes']['informative'].extend(positive_seqs)
                    sample_data[sample]['genotypes']['informative'].extend(partial_pos_seqs)
                if len(positive_seqs) == 1 and len(partial_pos_seqs) == 0:
                    sample_data[sample]['genotypes']['unique'].extend(positive_seqs)
                elif len(positive_seqs) == 0 and len(partial_pos_seqs) == 1:
                    sample_data[sample]['genotypes']['unique'].extend(partial_pos_seqs)


                #exclude genotypes which must have the negative kmer state
                if ratio == 1:
                    sample_data[sample]['genotypes']['exclude'].extend(list(set(genotypes) - set(positive_seqs + partial_pos_seqs)))


            elif sum(sample_data[sample]['counts'][target].values()) > detection_limit:
                sample_data[sample]['genotypes']['exclude'].extend(positive_seqs)

            sample_data[sample]['genotypes']['include'] = list(set(sample_data[sample]['genotypes']['include']))
            sample_data[sample]['genotypes']['exclude'] = list(set(sample_data[sample]['genotypes']['exclude']))
            sample_data[sample]['genotypes']['unique'] = list(set(sample_data[sample]['genotypes']['unique']))

    for sample in sample_data:
        include = sample_data[sample]['genotypes']['include']
        exclude = sample_data[sample]['genotypes']['exclude']
        informative = sample_data[sample]['genotypes']['informative']
        unique = sample_data[sample]['genotypes']['unique']

        if len(include) > 0:
            sample_data[sample]['genotypes']['candidates'] = list(set(include) - set(exclude))
            sample_data[sample]['genotypes']['candidates'].sort()
            sample_data[sample]['genotypes']['candidates'] = list(set(sample_data[sample]['genotypes']['candidates'] ) & set(informative))
            sample_data[sample]['genotypes']['candidates'].sort()


    return sample_data

def get_type_specific_targets(scheme_df):
    kmers = {}
    for row in scheme_df.itertuples():
        target = str(row.key)
        if isinstance(row.positive_seqs,float):
            positive_seqs = []
        else:
            positive_seqs = str(row.positive_seqs).split(',')

        if isinstance(row.partial_positive_seqs,float):
            partial_positive_seqs = []
        else:
            partial_positive_seqs = str(row.partial_positive_seqs).split(',')


        if len(positive_seqs) == 1 and len(partial_positive_seqs) == 0:
            kmers[target] = positive_seqs[0]
        if len(positive_seqs) == 0 and len(partial_positive_seqs) == 1:
            kmers[target] = partial_positive_seqs[0]

    return kmers

def get_pos_kmers(kmer_ratios,min_cov_frac=0.05):
    mixed_targets = []
    for target in kmer_ratios:
        ratio = kmer_ratios[target]
        if ratio >= min_cov_frac :
            mixed_targets.append(target)
    return mixed_targets

def type_kmer_mapping(scheme_df):
    kmers = {}
    for row in scheme_df.itertuples():
        target = str(row.key)
        if isinstance(row.positive_seqs,float):
            positive_seqs = []
        else:
            positive_seqs = str(row.positive_seqs).split(',')

        if isinstance(row.partial_positive_seqs,float) :
            partial_positive_seqs = []
        else:
            partial_positive_seqs = str(row.partial_positive_seqs).split(',')


        kmers[target] = positive_seqs + partial_positive_seqs

    return kmers


def type_occamization(sample_data,scheme_df,min_cov_frac=0.05,min_cov=20):
    type_kmer_mappings = type_kmer_mapping(scheme_df)
    type_specific_kmers = get_type_specific_targets(scheme_df)

    genotypes_with_specific_kmers = []
    for target in type_specific_kmers:
        genotypes_with_specific_kmers.append(type_specific_kmers[target])

    for genotype in ['B.1.1.51', 'B.1.1.51', 'B.1.1.100', 'B.1.292', 'A.28', 'A.18', 'B.1.1.7']:
        if genotype not in type_specific_kmers:
            print("{} no specific kmer".format(genotype))
        else:
            print("{} found specific kmer".format(genotype))

    for sample in sample_data:
        if not 'genotypes' in sample_data[sample]:
            continue

        type_data = sample_data[sample]['genotypes']
        candidates = type_data['candidates']
        for genotype in ['B.1.1.51','B.1.1.51','B.1.1.100','B.1.292','A.28','A.18','B.1.1.7']:
            if genotype not in candidates:
                print("{} not found".format(genotype))
            else:
                print("{} found".format(genotype))

        #skip samples where no genotype is compatible with the kmer data or a single candidate is available
        if len(candidates) <= 1:
            continue

        exclude = type_data['exclude']
        unique = type_data['unique']

        informative = type_data['informative']
        ratios = sample_data[sample]['ratios']
        counts = sample_data[sample]['counts']
        pos_sites = get_pos_kmers(ratios)
        simplified_type_set = [] + unique

        for genotype in ['B.1.1.51','B.1.1.100','B.1.292','A.28','A.18','B.1.1.7']:
            if genotype not in exclude:
                print("{} not excluded".format(genotype))
            else:
                print("{} excluded".format(genotype))

            if genotype not in informative:
                print("{} not informative".format(genotype))
            else:
                print("{} informative".format(genotype))

        genotype_info = {}

        for target in pos_sites:

            if not target in type_kmer_mappings:
                logging.warning("Cannot find target kmer {} in type mapping data".format(target))
                continue

            if ratios[target] < min_cov_frac or counts[target]['positive'] < min_cov:
                continue
            genotypes = type_kmer_mappings[target]
            compatible_genotypes = list(set(genotypes) - set(exclude))
            if len(compatible_genotypes) == 1 :
                print("{}\t{}".format(compatible_genotypes,target))
                unique.append(compatible_genotypes[0])
                simplified_type_set.append(compatible_genotypes[0])

            intersect = list(set(compatible_genotypes) & set(simplified_type_set ))

            if len(intersect) > 0:
                continue

            compatible_genotypes = list(set(compatible_genotypes) & set(informative ))

            for genotype in compatible_genotypes:
                simplified_type_set.append(genotype)


        genotype_unique_counts = {k: v for k, v in sorted(Counter(unique).items(), key=lambda item: item[1],reverse=True)}
        print(genotype_unique_counts)
        genotype_kmer_count = {}
        simplified_type_set = list(set(simplified_type_set))

        candidate_data = sample_data[sample]['genotypes']['candidate_data']

        for genotype in genotype_unique_counts:
            targets = candidate_data[genotype]['targets']
            print("{}\t{}".format(genotype,targets))
            if len(targets) == 0:
                continue
            for g2 in candidate_data :
                if g2  == genotype :
                    continue
                print(g2)
                candidate_data[g2]['targets'] = list(set(candidate_data[g2]['targets']) - set(targets))
        filtered = {}
        for genotype in candidate_data:
            if len(candidate_data[genotype]['targets']) > 0 and genotype in simplified_type_set:
                filtered[genotype] = candidate_data[genotype]
        candidate_data = filtered
        del(filtered)

        for genotype in candidate_data:
            genotype_kmer_count[genotype] = len(candidate_data[genotype]['targets'])
        genotype_target_counts = {k: v for k, v in sorted(genotype_kmer_count.items(), key=lambda item: item[1],reverse=True)}

        for genotype in genotype_target_counts:
            targets = candidate_data[genotype]['targets']
            print("{}\t{}".format(genotype,targets))
            if len(targets) == 0:
                continue
            for g2 in candidate_data :
                if g2  == genotype :
                    continue
                candidate_data[g2]['targets'] = list(set(candidate_data[g2]['targets']) - set(targets))
        filtered = {}
        for genotype in candidate_data:
            if len(candidate_data[genotype]['targets']) > 0 and genotype in simplified_type_set:
                filtered[genotype] = candidate_data[genotype]
        candidate_data = filtered



        sample_data[sample]['genotypes']['candidate_data'] = filtered
        sample_data[sample]['genotypes']['unique'] = list(set(unique) - set(list(candidate_data.keys())))
        sample_data[sample]['genotypes']['candidates'] = list(candidate_data.keys())



    return sample_data


def calc_type_coverage(sample_data,scheme_df,min_cov=20,min_cov_frac=0.05):
    genotypes_to_kmers = get_scheme_genotypes(scheme_df)

    for sample in sample_data:
        sample_data[sample]['genotypes']['candidate_data'] = {}
        if not 'genotypes' in sample_data[sample]:
            continue
        if not 'counts' in sample_data[sample]:
            continue
        if not 'ratios' in sample_data[sample]:
            continue

        type_data = sample_data[sample]['genotypes']
        candidates = type_data['candidates']
        frequencies = {}
        for genotype in candidates:
            if not genotype in genotypes_to_kmers:
                continue
            if not genotype in frequencies:
                frequencies[genotype] = {'data':{}}
            targets = genotypes_to_kmers[genotype]

            counts = []
            ratios = []
            for target in targets:
                if target not in sample_data[sample]['counts']:
                    continue
                count = sample_data[sample]['counts'][target]['positive']
                ratio = sample_data[sample]['ratios'][target]

                if count > min_cov:
                    counts.append(count)
                    if ratio > min_cov_frac:
                        ratios.append(ratio)
                        frequencies[genotype]['data'][target] = {'count':count,'ratios':ratio}


            sample_data[sample]['genotypes']['candidate_data'][genotype] ={'num_targets':0,'targets':list(frequencies[genotype]['data'].keys()),
                                                    'average_target_freq':0,
                                                    'target_freq_stdev':0,
                                                    'average_ratio':0,
                                                    'ratio_stdev':0}
            sample_data[sample]['genotypes']['candidate_data'][genotype]['num_targets'] = len(counts)
            if len(counts) == 0:
                continue
            total_freq = sum(counts)
            if total_freq > 0:
                sample_data[sample]['genotypes']['candidate_data'][genotype]['average_target_freq'] = total_freq / len(counts)
                sample_data[sample]['genotypes']['candidate_data'][genotype]['target_freq_stdev'] = statistics.pstdev(counts)
                sample_data[sample]['genotypes']['candidate_data'][genotype]['average_ratio'] = sum(ratios) / len(ratios)
                sample_data[sample]['genotypes']['candidate_data'][genotype]['ratio_stdev'] = statistics.pstdev(ratios)


    return sample_data



def get_scheme_genotypes(scheme_df):
    genotypes  = {}
    for row in scheme_df.itertuples():
        target = str(row.key)
        if isinstance(row.positive_seqs,float):
            positive_seqs = []
        else:
            positive_seqs = str(row.positive_seqs).split(',')

        if isinstance(row.partial_positive_seqs,float) :
            partial_positive_seqs = []
        else:
            partial_positive_seqs = str(row.partial_positive_seqs).split(',')

        to_add = positive_seqs + partial_positive_seqs
        for genotype in to_add:
            if not genotype in genotypes:
                genotypes[genotype] = []
            genotypes[genotype].append(target)
    return genotypes




def main():
    cmd_args = parse_args()
    logger = init_console_logger(2)
    is_args_ok = validate_args(cmd_args,logger)
    if not is_args_ok:
        logger.Error("One or more command line arguments has an issue, please check the log messages and try again")
        sys.exit()

    #input parameters
    prefix = cmd_args.prefix
    min_cov = cmd_args.min_cov
    min_cov_frac = cmd_args.min_cov_frac
    detection_limit = 100
    max_missing_sites = cmd_args.max_missing_sites
    max_mixed_sites = cmd_args.max_mixed_sites
    scheme = cmd_args.scheme
    mode = cmd_args.mode
    type = cmd_args.mode
    R1 = cmd_args.R1
    R2 = cmd_args.R2
    SE = cmd_args.se
    data_dir = cmd_args.data_dir
    outdir = cmd_args.outdir
    nthreads = cmd_args.n_threads
    strain_profiles = cmd_args.strain_profiles
    genotype_profiles = cmd_args.genotype_profiles
    primers = cmd_args.primers
    no_template_control = cmd_args.no_template_control
    positive_control = cmd_args.positive_control
    genotype_dist_cutoff = cmd_args.genotype_dist_cutoff
    genotype_dist_cutoff = cmd_args.genotype_dist_cutoff

    # initialize analysis directory
    if not os.path.isdir(outdir):
        logger.info("Creating analysis results directory {}".format(outdir))
        os.mkdir(outdir, 0o755)
    else:
        logger.info("Results directory {} already exits, will overwrite any results files here".format(outdir))

    logger.info("Reading kmer scheme from {}".format(scheme))
    if scheme in TYPING_SCHEMES:
        scheme = TYPING_SCHEMES[scheme]

    scheme_df = read_tsv(scheme)
    logger.info("Found {} kmers in {}".format(len(scheme_df),scheme))
    biohansel_fasta_file = os.path.join(outdir, "{}.biohansel.fasta".format(prefix))
    logger.info("Writing biohansel compatible kmer scheme to {}".format(biohansel_fasta_file))

    scheme_to_biohansel_fasta(scheme_df,biohansel_fasta_file )
    scheme_kmer_target_info = init_kmer_targets(scheme_df)
    scheme_kmer_target_keys = list(scheme_kmer_target_info.keys())
    scheme_kmer_target_keys.sort()
    scheme_kmer_target_keys = [str(i) for i in scheme_kmer_target_keys]
    #scheme_kmer_names = generate_biohansel_kmer_names(scheme_kmer_target_keys)

    #Identify kmers which are present in the no template control for masking purposes later
    no_template_kmers = {}
    if no_template_control is not None:
        kmer_file = os.path.join(outdir,"{}.no_template.bh.kmer.txt")
        summary_file = os.path.join(outdir,"{}.no_template.bh.summary.txt")
        simple_file = os.path.join(outdir,"{}.no_template.bh.simple.txt")
        logger.info("Identifying kmers which are found in no template control {}".format(no_template_control))
        bio_hansel.run_biohansel_single(biohansel_fasta_file, no_template_control, kmer_file, summary_file, simple_file,
                                 min_cov, nthreads)
        no_template_df = read_tsv(kmer_file)
        no_template_df = no_template_df[no_template_df['freq'] >= min_cov]
        no_template_kmers = no_template_df[['kmername','freq']].to_dict()


    #Identify kmers which are present in the primer sets for masking purposes later
    primer_kmers = {}
    if primers is not None:
        if primers in PRIMER_SCHEMES:
            primers = PRIMER_SCHEMES[primers]
        logger.info("Writing primers to fasta file".format(primers))
        primer_df = read_tsv(primers)
        primer_file = os.path.join(outdir,"primers.fasta")
        fh = open(primer_file,'w')
        for row in primer_df.itertuples():
            fh.write(">{}\n{}\n".format(row.name,row.seq))
        fh.close()

        kmer_file = os.path.join(outdir,"primers.bh.kmer.txt")
        summary_file = os.path.join(outdir,"primers.bh.summary.txt")
        simple_file = os.path.join(outdir,"primers.bh.simple.txt")
        logger.info("Identifying kmers which are found in selected primer set {}".format(primers))
        bio_hansel.run_biohansel_single(biohansel_fasta_file, primer_file, kmer_file, summary_file, simple_file,
                                 min_cov, nthreads)
        primer_df = read_tsv(kmer_file)
        primer_kmers = {}
        for kmer in primer_df['kmername'].to_list():
            primer_kmers[kmer] = min_cov

    #Identify kmers identified in positive control
    positive_control_kmers = {}
    if positive_control is not None:
        kmer_file = os.path.join(outdir, "positive_control.bh.kmer.txt")
        summary_file = os.path.join(outdir, "positive_control.bh.summary.txt")
        simple_file = os.path.join(outdir, "positive_control.bh.simple.txt")
        logger.info("Identifying kmers which are found in positive control {}".format(positive_control))
        bio_hansel.run_biohansel_single(biohansel_fasta_file, positive_control, kmer_file, summary_file, simple_file,
                                        min_cov, nthreads)
        positive_control_df = read_tsv(kmer_file)
        positive_control_df = positive_control_df[positive_control_df['freq'] >= min_cov]
        positive_control_kmers = positive_control_df[['kmername', 'freq']].to_dict()

    #Run kmer detection on the samples
    if data_dir is not None:
        kmer_file = os.path.join(outdir,"{}.sample.bh.kmer.txt".format(prefix))
        summary_file = os.path.join(outdir,"{}.sample.bh.summary.txt".format(prefix))
        simple_file = os.path.join(outdir,"{}.sample.bh.simple.txt".format(prefix))
        logger.info("Identifying kmers which are found in input {}".format(data_dir))
        bio_hansel.run_biohansel_directory(biohansel_fasta_file, data_dir, kmer_file, summary_file, simple_file,
                                 min_cov, nthreads)
    elif SE is not None:
        kmer_file = os.path.join(outdir,"{}.sample.bh.kmer.txt".format(prefix))
        summary_file = os.path.join(outdir,"{}.sample.bh.summary.txt".format(prefix))
        simple_file = os.path.join(outdir,"{}.sample.bh.simple.txt".format(prefix))
        logger.info("Identifying kmers which are found in input {}".format(SE))
        bio_hansel.run_biohansel_single(biohansel_fasta_file, SE, kmer_file, summary_file, simple_file,
                                 min_cov, nthreads)
    else:
        kmer_file = os.path.join(outdir,"{}.sample.bh.kmer.txt".format(prefix))
        summary_file = os.path.join(outdir,"{}.sample.bh.summary.txt".format(prefix))
        simple_file = os.path.join(outdir,"{}.sample.bh.simple.txt".format(prefix))
        logger.info("Identifying kmers which are found in input {} {}".format(R1,R2))
        bio_hansel.run_biohansel_paired(biohansel_fasta_file, R1, R2, kmer_file, summary_file, simple_file,
                                 min_cov, nthreads)

    sample_kmer_biohansel_df = read_tsv(kmer_file)
    sample_kmer_data = calc_kmer_ratio(sample_kmer_biohansel_df,scheme_kmer_target_keys,min_cov)
    sample_kmer_data = calc_mixed_sites(sample_kmer_data,min_cov_frac)



    #subtract contamination kmers in the no template control
    if len(no_template_kmers) > 0:
        sample_kmer_data = subtract_contamination_kmers(sample_kmer_data, no_template_kmers, min_cov_frac)

    #subtract contamination kmers in the primers
    if len(primer_kmers) > 0:
        sample_kmer_data = subtract_contamination_kmers(sample_kmer_data, primer_kmers, min_cov_frac)


    #Get information about each of the targets involved in mixed events
    target_info = get_mixed_kmer_targets(sample_kmer_data, min_cov_frac)

    #Get list of strains which are compatible with the kmer information
    sample_kmer_data =identify_compatible_types(scheme_df, sample_kmer_data, min_cov_frac, detection_limit=100)

    sample_kmer_data = calc_type_coverage(sample_kmer_data, scheme_df)
    sample_kmer_data = type_occamization(sample_kmer_data,scheme_df)
    #sample_kmer_data = calc_type_coverage(sample_kmer_data, scheme_df)


    for sample in sample_kmer_data:

        total_ratio = 0
        sample_kmer_data[sample]['genotypes']['candidate_data']
        for genotype in sample_kmer_data[sample]['genotypes']['candidate_data']:
            data = sample_kmer_data[sample]['genotypes']['candidate_data'][genotype]

            #Filter out strains with low confidence
            if data['average_ratio'] <= min_cov_frac or data['average_target_freq'] <= min_cov:
                continue
            data['targets'].sort()
            total_ratio += data['average_ratio']
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(sample,
                                      genotype,
                                      data['num_targets'],
                                      data['targets'],
                                      data['average_target_freq'],
                                      data['target_freq_stdev'],
                                      data['average_ratio'],
                                      data['ratio_stdev']))
        unknown_fraq = 1 - total_ratio
        if unknown_fraq >= min_cov_frac:
            print("{}\tunknown\t0\t0\t0\t{}\t{}".format(sample,unknown_fraq,0))


    # call main function
if __name__ == '__main__':
    main()