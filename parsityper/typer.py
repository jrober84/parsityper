#!/usr/bin/python
import time
from argparse import (ArgumentParser, FileType)
import logging, os, sys, re, collections, operator, math, shutil, datetime
from collections import Counter
import pandas as pd
from parsityper.helpers import validate_args, init_console_logger, read_tsv, scheme_to_biohansel_fasta, \
    filter_biohansel_kmer_df, process_biohansel_kmer, \
    init_kmer_targets, generate_biohansel_kmer_names, get_scheme_template, generate_random_phrase, calc_md5, \
    get_kmer_groups, get_kmer_group_mapping, get_sequence_files, read_genotype_profiles, dist_compatible_profiles, \
    generate_target_presence_table, nan_compatible_kmer_pairwise_distmatrix
import copy
from parsityper.bio_hansel import bio_hansel
import statistics
from parsityper.constants import PRIMER_SCHEMES, TYPING_SCHEMES
from parsityper.helpers import profile_pairwise_distmatrix
from parsityper.visualizations import dendrogram_visualization, plot_mds, generate_sample_coverage_plot


def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Parsimony based sample Kmer-analysis')
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
    parser.add_argument('--scheme', type=str, required=True,
                        help='TSV formated kmer scheme', default='')
    parser.add_argument('--primers', type=str, required=False,
                        help='TSV formated primer file', default=None)
    parser.add_argument('--outdir', type=str, required=True,
                        help='output directory')
    parser.add_argument('--prefix', type=str, required=False,
                        help='output file prefix', default='SARS-COV-2_kmer_analysis')
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
    parser.add_argument('--genotype_profiles', type=str, required=False,
                        help='Kmer profile scheme to compare against')
    parser.add_argument('--no_template_control', type=str, required=False,
                        help='Fasta/Fastq formatted no template control data')
    parser.add_argument('--positive_control', type=str, required=False,
                        help='Fasta/Fastq formatted positive control data')
    parser.add_argument('--list', type=str, required=False,
                        help='list in-built primer and typing schemes')
    return parser.parse_args()


def generate_sample_kmer_profile(sample_kmer_data, scheme_kmer_target_keys, min_cov=20):
    '''
    Accepts biohansel kmer pandas dataframe and sums the frequency of positive and negative kmers by sample if they meet the minimum threshold
    :param biohansel_df: pandas biohansel kmer df
    :param scheme_kmer_target_keys: list of kmer targets
    :param min_cov: integer of minimum frequency
    :return: dict of samples and counts of positive and negative kmers
    '''
    samples = {}

    for sample in sample_kmer_data:
        if not sample in samples:
            samples[sample] = {}
            value = {'positive': 0, 'negative': 0}
            samples[sample]['counts'] = copy.deepcopy(get_scheme_template(scheme_kmer_target_keys, value))
        for target in sample_kmer_data[sample]:
            if not target in samples[sample]['counts']:
                continue
            for seq in sample_kmer_data[sample][target]:
                is_pos_kmer = sample_kmer_data[sample][target][seq]['is_pos_kmer']
                freq = sample_kmer_data[sample][target][seq]['freq']

                if freq >= min_cov:
                    if is_pos_kmer:
                        samples[sample]['counts'][target]['positive'] += freq
                    else:
                        samples[sample]['counts'][target]['negative'] += freq

    return samples


def calc_kmer_ratio(sample_kmer_data, scheme_kmer_target_keys, min_cov=20):
    '''
    Accepts biohansel kmer pandas dataframe and sums the frequency of positive and negative kmers by sample if they meet the minimum threshold
    and then determines the ratio of the positive and negative kmers
    :param biohansel_df: pandas biohansel kmer df
    :param scheme_kmer_target_keys: list of kmer targets
    :param min_cov: integer of minimum frequency
    :return: dict of samples and counts of positive and negative kmers
    '''
    samples = generate_sample_kmer_profile(sample_kmer_data, scheme_kmer_target_keys, min_cov)

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


def calc_mixed_sites(sample_data, min_cov_frac):
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
            sample_data[sample]['statistics']['num_kmers_targets_present'] += 1
            if ratio == 1 or ratio == 0:
                continue

            if ratio >= min_cov_frac or ratio <= (1 - min_cov_frac):
                sample_data[sample]['statistics']['num_mixed_kmers_targets'] += 1
                sample_data[sample]['statistics']['mixed_kmers_targets'].append(target)
                mixed_ratios.append(ratio)
        if len(total_freqs) > 0:
            sample_data[sample]['statistics']['average_target_coverage'] = sum(total_freqs) / len(total_freqs)
            sample_data[sample]['statistics']['stdev_target_coverage'] = statistics.pstdev(total_freqs)
        else:
            sample_data[sample]['statistics']['average_target_coverage'] = 0
            sample_data[sample]['statistics']['stdev_target_coverage'] = 0

        if len(mixed_ratios) > 0:
            sample_data[sample]['statistics']['average_mixed_target_ratio'] = sum(mixed_ratios) / len(mixed_ratios)
            sample_data[sample]['statistics']['stdev_mixed_target_ratio'] = statistics.pstdev(mixed_ratios)
        else:
            sample_data[sample]['statistics']['average_mixed_target_ratio'] = 0
            sample_data[sample]['statistics']['stdev_mixed_target_ratio'] = 0

    return sample_data


def get_mixed_kmer_targets(sample_data, min_cov_frac):
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
                    'count_samples_present': 0,
                    'count_samples_absent': 0,
                    'count_mixed_samples': 0,
                    'mixed_samples_ids': []
                }
            ratio = sample_data[sample]['ratios'][target]
            if ratio == -1:
                targets[target]['count_samples_absent'] += 1
            else:
                targets[target]['count_samples_present'] += 1
                if (ratio != 0 and ratio != 1) and (ratio >= min_cov_frac or ratio <= (1 - min_cov_frac)):
                    targets[target]['count_mixed_samples'] += 1
                    targets[target]['mixed_samples_ids'].append(sample)
    return targets


def subtract_contamination_kmers(sample_data, kmers, min_cov_frac):
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
            target = kmer.replace('negative', '').split('-')[0]
            frequency = kmers[kmer]
            # The level of contamination in this sample is below what would be called at a sub-consensus level so no need to change
            # anything
            if target not in sample_data[sample]['ratios']:
                continue
            if sample_data[sample]['ratios'][target] < min_cov_frac or sample_data[sample]['ratios'][
                target] > 1 - min_cov_frac:
                continue
            sample_data[sample]['counts'][target][id] -= frequency
            if sample_data[sample]['counts'][target][id] < 0:
                sample_data[sample]['counts'][target][id] = 0
            total = sum(sample_data[sample]['counts'][target].values())
            if total == 0:
                sample_data[sample]['ratios'][target] = -1
            else:
                sample_data[sample]['ratios'][target] = sample_data[sample]['counts'][target]['positive'] / total
    return sample_data


def identify_compatible_types(scheme_df, sample_data, min_cov_frac, detection_limit=100):
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

    for sample in sample_data:
        sample_data[sample]['genotypes'] = {
                'include': [],
                'informative': [],
                'exclude': [],
                'unique': [],
                'candidates': [],
                'candidate_data': {}

            }

    for row in scheme_df.itertuples():

        if isinstance(row.positive_seqs,float):
            positive_seqs = []
        else:
            positive_seqs = str(row.positive_seqs).split(',')
        if isinstance(row.partial_positive_seqs,float):
            partial_pos_seqs = []
        else:
            partial_pos_seqs = str(row.partial_positive_seqs).split(',')

        target = str(row.key)

        for sample in sample_data:

            if target not in sample_data[sample]['counts']:
                continue

            ratio = sample_data[sample]['ratios'][target]
            if ratio == -1:
                continue
            if ratio >= min_cov_frac:
                sample_data[sample]['genotypes']['include'].extend(positive_seqs)
                if len(partial_pos_seqs) < len(genotypes):
                    sample_data[sample]['genotypes']['include'].extend(partial_pos_seqs)
                if len(positive_seqs) < len(genotypes) * 0.01 and len(partial_pos_seqs) < len(genotypes) * 0.01:
                    sample_data[sample]['genotypes']['informative'].extend(positive_seqs)
                    sample_data[sample]['genotypes']['informative'].extend(partial_pos_seqs)
                if len(positive_seqs) == 1 and len(partial_pos_seqs) == 0:
                    sample_data[sample]['genotypes']['unique'].extend(positive_seqs)
                elif len(positive_seqs) == 0 and len(partial_pos_seqs) == 1:
                    sample_data[sample]['genotypes']['unique'].extend(partial_pos_seqs)
                # exclude genotypes which must have the negative kmer state
                # This needs more thought as it has caused some issues
                if ratio >= 1 - min_cov_frac and sum(sample_data[sample]['counts'][target].values()) >= detection_limit:
                    if len(positive_seqs) > 0:
                        e_list = list(set(genotypes) - set(positive_seqs + partial_pos_seqs))
                        sample_data[sample]['genotypes']['exclude'].extend(e_list)

            elif sum(sample_data[sample]['counts'][target].values()) >= detection_limit:
               sample_data[sample]['genotypes']['exclude'].extend(positive_seqs)

            sample_data[sample]['genotypes']['include'] = list(
                set(sample_data[sample]['genotypes']['include']) - set(['nan']))
            if len(sample_data[sample]['genotypes']['include']) == 0:
                sample_data[sample]['genotypes']['include'] = list(
                    set(genotypes) - set(sample_data[sample]['genotypes']['exclude']))
            sample_data[sample]['genotypes']['exclude'] = list(
                set(sample_data[sample]['genotypes']['exclude']) - set(['nan']))
            sample_data[sample]['genotypes']['unique'] = list(
                set(sample_data[sample]['genotypes']['unique']) - set(['nan']))
            sample_data[sample]['genotypes']['informative'] = list(
                set(sample_data[sample]['genotypes']['informative']) - set(['nan']))

    for sample in sample_data:
        include = sample_data[sample]['genotypes']['include']
        exclude = sample_data[sample]['genotypes']['exclude']


        if len(include) > 0:
            sample_data[sample]['genotypes']['candidates'] = list(set(include) - set(exclude))
            sample_data[sample]['genotypes']['candidates'].sort()
            sample_data[sample]['genotypes']['candidates'] = list(
                set(sample_data[sample]['genotypes']['candidates']) - set('nan'))

    return sample_data


def get_type_specific_targets(scheme_df):
    '''

    :param scheme_df:
    :type scheme_df:
    :return:
    :rtype:
    '''
    kmers = {}
    for row in scheme_df.itertuples():
        target = str(row.key)
        if isinstance(row.positive_seqs, float):
            positive_seqs = []
        else:
            positive_seqs = str(row.positive_seqs).split(',')

        if isinstance(row.partial_positive_seqs, float):
            partial_positive_seqs = []
        else:
            partial_positive_seqs = str(row.partial_positive_seqs).split(',')

        if len(positive_seqs) == 1 and len(partial_positive_seqs) == 0:
            kmers[target] = positive_seqs[0]
        if len(positive_seqs) == 0 and len(partial_positive_seqs) == 1:
            kmers[target] = partial_positive_seqs[0]

    return kmers


def get_pos_kmers(kmer_ratios, min_cov_frac=0.05):
    '''

    :param kmer_ratios:
    :type kmer_ratios:
    :param min_cov_frac:
    :type min_cov_frac:
    :return:
    :rtype:
    '''
    mixed_targets = []
    for target in kmer_ratios:
        ratio = kmer_ratios[target]
        if ratio >= min_cov_frac:
            mixed_targets.append(target)
    return mixed_targets


def type_kmer_mapping(scheme_df):
    '''

    :param scheme_df:
    :type scheme_df:
    :return:
    :rtype:
    '''
    kmers = {}
    for row in scheme_df.itertuples():
        target = str(row.key)
        if isinstance(row.positive_seqs, float):
            positive_seqs = []
        else:
            positive_seqs = str(row.positive_seqs).split(',')

        if isinstance(row.partial_positive_seqs, float):
            partial_positive_seqs = []
        else:
            partial_positive_seqs = str(row.partial_positive_seqs).split(',')

        kmers[target] = positive_seqs + partial_positive_seqs

    return kmers

def assign_targets(candidate_data,genotype_kmer_count):
    genotype_target_counts = {k: v for k, v in
                              sorted(genotype_kmer_count.items(), key=lambda item: item[1], reverse=True)}
    delta = False

    for genotype in genotype_target_counts:
        targets = candidate_data[genotype]['targets']
        for g2 in candidate_data:
            if g2 == genotype:
                continue
            new = list(set(candidate_data[g2]['targets']) - set(targets))
            old = list(set(candidate_data[g2]['targets']) )
            candidate_data[g2]['targets'] = new
            candidate_data[g2]['num_targets'] = len(new)
            if len(new) != len(old):
                delta = True

        if delta:
            break
    return (delta,candidate_data)

def get_genotype_target_associations(scheme_df):
    '''
    :param scheme_df:
    :type scheme_df:
    :return:
    :rtype:
    '''
    genotypes = {}
    for row in scheme_df.itertuples():
        target = str(row.key)
        if isinstance(row.positive_seqs, float):
            positive_seqs = []
        else:
            positive_seqs = str(row.positive_seqs).split(',')

        if isinstance(row.partial_positive_seqs, float):
            partial_positive_seqs = []
        else:
            partial_positive_seqs = str(row.partial_positive_seqs).split(',')
        for g in positive_seqs:
            if not g in genotypes:
                genotypes[g] = {'positive':[],'partial':[]}
            genotypes[g]['positive'].append(target)
        for g in partial_positive_seqs:
            if not g in genotypes:
                genotypes[g] = {'positive':[],'partial':[]}
            genotypes[g]['partial'].append(target)
    return genotypes


def type_occamization(sample_data, scheme_df, min_cov_frac=0.05, min_cov=20):
    type_kmer_mappings = type_kmer_mapping(scheme_df)
    type_specific_kmers = get_type_specific_targets(scheme_df)
    genotype_kmer_assoc = get_genotype_target_associations(scheme_df)

    genotypes_with_specific_kmers = []
    for target in type_specific_kmers:
        genotypes_with_specific_kmers.append(type_specific_kmers[target])

    for sample in sample_data:
        if not 'genotypes' in sample_data[sample]:
            logging.warning("Could not find genoytpe data for sample {}...skip".format(sample))
            continue

        type_data = sample_data[sample]['genotypes']
        candidates = type_data['candidates']

        # skip samples where no genotype is compatible with the kmer data or a single candidate is available
        if len(candidates) <= 1:
            continue

        exclude = type_data['exclude']
        unique = type_data['unique']

        informative = type_data['informative']
        ratios = sample_data[sample]['ratios']
        counts = sample_data[sample]['counts']
        pos_sites = get_pos_kmers(ratios)

        if len(pos_sites) == 0:
            candidate_data = sample_data[sample]['genotypes']['candidate_data']
            for genotype in candidates:
                num_targets = 1 + len(candidate_data[genotype]['targets'])
                num_positives = 1 + len(genotype_kmer_assoc[genotype]['positive'])
            continue

        simplified_type_set = [] + unique

        for target in pos_sites:
            if not target in type_kmer_mappings:
                logging.warning("Cannot find target kmer {} in type mapping data".format(target))
                continue

            if ratios[target] < min_cov_frac or counts[target]['positive'] < min_cov:
                logging.info("{} target kmer {} does not meet coverage {} and ratio {} requirement".format(sample,
                                                                                                           target,
                                                                                                           counts[
                                                                                                               target][
                                                                                                               'positive'],
                                                                                                           ratios[
                                                                                                               target]))
                continue
            genotypes = type_kmer_mappings[target]
            compatible_genotypes = list(set(genotypes) - set(exclude))
            if len(compatible_genotypes) == 1:
                unique.append(compatible_genotypes[0])
                simplified_type_set.append(compatible_genotypes[0])

            intersect = list(set(compatible_genotypes) & set(simplified_type_set))

            if len(intersect) > 0:
                continue

            # compatible_genotypes = list(set(compatible_genotypes) & set(informative ))

            for genotype in compatible_genotypes:
                simplified_type_set.append(genotype)

        genotype_unique_counts = {k: v for k, v in
                                  sorted(Counter(unique).items(), key=lambda item: item[1], reverse=True)}

        genotype_kmer_count = {}
        simplified_type_set = list(set(simplified_type_set))
        candidate_data = sample_data[sample]['genotypes']['candidate_data']

        for i in genotype_unique_counts:
            genotype = genotype_unique_counts[i]
            if genotype not in candidate_data:
                continue
            targets = candidate_data[genotype]['targets']
            if len(targets) == 0:
                continue
            for k in range(i + 1, len(genotype_unique_counts)):
                g2 = genotype_unique_counts[k]
                candidate_data[g2]['targets'] = list(set(candidate_data[g2]['targets']) - set(targets))

        filtered = {}
        for genotype in candidate_data:
            if len(candidate_data[genotype]['targets']) > 0 and genotype in simplified_type_set:
                filtered[genotype] = candidate_data[genotype]
        candidate_data = filtered
        del (filtered)

        found_pos_kmer_targets = []
        for genotype in candidate_data:
            genotype_kmer_count[genotype] = len(candidate_data[genotype]['targets'])
            found_pos_kmer_targets.extend(candidate_data[genotype]['targets'])
        genotype_target_counts = {k: v for k, v in
                                  sorted(genotype_kmer_count.items(), key=lambda item: item[1], reverse=True)}

        (status,candidate_data) = assign_targets(candidate_data, genotype_kmer_count)

        while status:
            for genotype in candidate_data:
                genotype_kmer_count[genotype] = len(candidate_data[genotype]['targets'])
            (status,candidate_data) = assign_targets(candidate_data, genotype_kmer_count)


        filtered = {}
        for genotype in candidate_data:
            if len(candidate_data[genotype]['targets']) > 0 and genotype in simplified_type_set:
                filtered[genotype] = candidate_data[genotype]
        candidate_data = filtered
        sample_data[sample]['genotypes']['candidate_data'] = filtered
        sample_data[sample]['genotypes']['unique'] = list(set(unique) - set(list(candidate_data.keys())))
        sample_data[sample]['genotypes']['candidates'] = list(candidate_data.keys())
    return sample_data


def get_count_positive_targets(counts, min_cov):
    count_pos = 0
    for target in counts:
        if counts[target]['positive'] >= min_cov:
            count_pos += 1
    return count_pos


def calc_type_coverage(sample_data, scheme_df, min_cov=20, min_cov_frac=0.05, recalc=False):
    genotypes_to_kmers = get_scheme_genotypes(scheme_df)
    for sample in sample_data:
        if not recalc:
            sample_data[sample]['genotypes']['candidate_data'] = {}
        if not 'genotypes' in sample_data[sample]:
            continue
        if not 'counts' in sample_data[sample]:
            continue
        if not 'ratios' in sample_data[sample]:
            continue
        count_pos = get_count_positive_targets(sample_data[sample]['counts'], min_cov)

        type_data = sample_data[sample]['genotypes']
        candidates = type_data['candidates']
        if 'nan' in candidates:
            candidates = list(set(candidates) - set(['nan']))
        if recalc:
            candidates = sample_data[sample]['genotypes']['candidates']
        frequencies = {}

        for genotype in candidates:
            if not genotype in genotypes_to_kmers:
                logging.warning("Warning not genotype to kmer mapping for type {}".format(genotype))
                continue
            if not genotype in frequencies:
                frequencies[genotype] = {'data': {}}
            targets = list(set(genotypes_to_kmers[genotype]))

            counts = []
            ratios = []

            # Accounts for genotypes which have no postitive targets
            count_targets = len(targets)

            if count_pos == 0:
                counts = sample_data[sample]['counts']
                ratios = sample_data[sample]['ratios']
                sample_counts = []
                sample_ratios = []
                for target in counts:
                    neg = counts[target]['negative']
                    pos = counts[target]['positive']
                    if neg >= min_cov or pos >= min_cov:
                        sample_counts.append(neg)
                        sample_ratios.append(ratios[target])

                sample_data[sample]['genotypes']['candidate_data'][genotype] = {'num_targets': 0, 'targets': list(
                    frequencies[genotype]['data'].keys()), 'average_target_freq': 0,
                                                                                'target_freq_stdev': 0,
                                                                                'average_ratio': 0,
                                                                                'ratio_stdev': 0}

                if len(counts) == 0:
                    continue
                total_freq = sum(sample_counts)
                if total_freq > 0:
                    sample_data[sample]['genotypes']['candidate_data'][genotype]['num_targets'] = len(sample_counts)
                    sample_data[sample]['genotypes']['candidate_data'][genotype][
                        'average_target_freq'] = total_freq / len(sample_counts)
                    sample_data[sample]['genotypes']['candidate_data'][genotype][
                        'target_freq_stdev'] = statistics.pstdev(sample_counts)
                    sample_data[sample]['genotypes']['candidate_data'][genotype]['average_ratio'] = sum(
                        sample_ratios) / len(sample_ratios)
                    sample_data[sample]['genotypes']['candidate_data'][genotype]['ratio_stdev'] = statistics.pstdev(
                        sample_ratios)

                continue

            for target in targets:
                if target not in sample_data[sample]['counts']:
                    logging.error("Warning target not in counts, data structure error {} - {}".format(genotype, target))
                    continue
                if recalc:
                    if target not in sample_data[sample]['genotypes']['candidate_data'][genotype]['targets']:
                        continue
                count = sample_data[sample]['counts'][target]['positive']
                ratio = sample_data[sample]['ratios'][target]

                if count >= min_cov:
                    if ratio >= min_cov_frac:
                        counts.append(count)
                        ratios.append(ratio)
                        frequencies[genotype]['data'][target] = {'count': count, 'ratios': ratio}

            sample_data[sample]['genotypes']['candidate_data'][genotype] = {'num_targets': 0, 'targets': list(
                frequencies[genotype]['data'].keys()),
                                                                            'average_target_freq': 0,
                                                                            'target_freq_stdev': 0,
                                                                            'average_ratio': 0,
                                                                            'ratio_stdev': 0}
            sample_data[sample]['genotypes']['candidate_data'][genotype]['num_targets'] = len(counts)

            if len(counts) == 0:
                continue
            total_freq = sum(counts)
            if total_freq > 0:
                sample_data[sample]['genotypes']['candidate_data'][genotype]['average_target_freq'] = total_freq / len(
                    counts)
                sample_data[sample]['genotypes']['candidate_data'][genotype]['target_freq_stdev'] = statistics.pstdev(
                    counts)
                sample_data[sample]['genotypes']['candidate_data'][genotype]['average_ratio'] = sum(ratios) / len(
                    ratios)
                sample_data[sample]['genotypes']['candidate_data'][genotype]['ratio_stdev'] = statistics.pstdev(ratios)

    return sample_data


def get_scheme_genotypes(scheme_df):
    '''

    :param scheme_df:
    :type scheme_df:
    :return:
    :rtype:
    '''
    genotypes = {}
    for row in scheme_df.itertuples():
        target = str(row.key)
        if isinstance(row.positive_seqs, float):
            positive_seqs = []
        else:
            positive_seqs = str(row.positive_seqs).split(',')

        if isinstance(row.partial_positive_seqs, float):
            partial_positive_seqs = []
        else:
            partial_positive_seqs = str(row.partial_positive_seqs).split(',')

        to_add = positive_seqs + partial_positive_seqs
        for genotype in to_add:
            if not genotype in genotypes:
                genotypes[genotype] = []
            genotypes[genotype].append(target)
    return genotypes


def write_sample_detailed_report(sample_report, outfile, scheme_kmer_target_info):
    '''

    :param sample_report:
    :type sample_report:
    :param outfile:
    :type outfile:
    :param scheme_kmer_target_info:
    :type scheme_kmer_target_info:
    :return:
    :rtype:
    '''
    report = ["\t".join(
        ['sample_id', 'genoytpe', 'num_targets', 'average_target_freq', 'target_freq_stdev', 'average_ratio',
         'ratio_stdev', 'assigned_positive_targets'])]
    for sample_id in sample_report:
        data = sample_report[sample_id]
        for genotype in data:
            targets = []
            if isinstance(data[genotype]['targets'], list):
                t = [int(x) for x in data[genotype]['targets']]
                t.sort()
            for target in data[genotype]['targets']:
                targets.append(scheme_kmer_target_info[target]['dna_name'])
            num_targets = len(data[genotype]['targets'])
            if num_targets == 0:
                num_targets = data[genotype]['num_targets']
            report.append("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(sample_id,
                                                                  genotype,
                                                                  num_targets,
                                                                  data[genotype]['average_target_freq'],
                                                                  data[genotype]['target_freq_stdev'],
                                                                  data[genotype]['average_ratio'],
                                                                  data[genotype]['ratio_stdev'],
                                                                  ', '.join(targets)))

    fh = open(outfile, 'w')
    fh.write("\n".join(report))


def get_detected_target_summary(sample_kmer_data, min_cov=20, min_cov_frac=0.05):
    '''

    :param sample_kmer_data:
    :type sample_kmer_data:
    :param min_cov:
    :type min_cov:
    :param min_cov_frac:
    :type min_cov_frac:
    :return:
    :rtype:
    '''
    targets = {}
    for sample_id in sample_kmer_data:
        if not sample_id in targets:
            targets[sample_id] = {'counts': [], 'positive': [], 'negative': [], 'targets': [], 'mixed': []}
        for target in sample_kmer_data[sample_id]['ratios']:
            ratio = sample_kmer_data[sample_id]['ratios'][target]
            count = sum(sample_kmer_data[sample_id]['counts'][target].values())
            if count >= min_cov:
                targets[sample_id]['counts'].append(count)
                targets[sample_id]['targets'].append(target)
                if ratio >= min_cov_frac:
                    targets[sample_id]['positive'].append(target)
                    if ratio < 1 - min_cov_frac and ratio > 0:
                        targets[sample_id]['mixed'].append(target)
                else:
                    targets[sample_id]['negative'].append(target)
    return targets


def write_sample_summary_report(sample_report, outfile, sample_kmer_data, scheme_kmer_target_info, sample_complexity,genes,
                                max_mixed_sites=10, min_cov=20, min_cov_frac=0.05,sample_type='single'):
    '''

    :param sample_report:
    :type sample_report:
    :param outfile:
    :type outfile:
    :param sample_kmer_data:
    :type sample_kmer_data:
    :param scheme_kmer_target_info:
    :type scheme_kmer_target_info:
    :param sample_complexity:
    :type sample_complexity:
    :param max_mixed_sites:
    :type max_mixed_sites:
    :param min_cov:
    :type min_cov:
    :param min_cov_frac:
    :type min_cov_frac:
    :return:
    :rtype:
    '''
    kmer_summary = get_detected_target_summary(sample_kmer_data, min_cov=min_cov, min_cov_frac=min_cov_frac)
    profile_names = {}
    header = ['sample_id','compatible_genotype(s)','sample_type','primary_genotype',
             'primary_genotype_transition_ratio','primary_genotype_transition_slope','num_targets_detected',
              'average_target_freq','target_freq_stdev','num_mixed_targets','kmer_profile_st','kmer_profile_md5',
              'kmer_profile_name','qc_sample_type','qc_max_missing_sites','qc_max_het_sites','qc_cov_stdev','qc_overall']
    for gene in genes:
        header.append("detected_mutations_{}_dna".format(gene))
    for gene in genes:
        header.append("detected_mutations_{}_aa".format(gene))
    report = ["\t".join(header)]

    for sample_id in sample_report:
        data = sample_report[sample_id]
        quality = sample_kmer_data[sample_id]['quality']
        mutations_dna = {}
        mutations_aa = {}
        positive_kmers = kmer_summary[sample_id]['positive']
        positive_kmers.sort()
        kmer_profile_st = positive_kmers
        for target in positive_kmers:

            dna_name = str(scheme_kmer_target_info[target]['dna_name'])
            aa_name = str(scheme_kmer_target_info[target]['aa_name'])
            gene = scheme_kmer_target_info[target]['gene']
            if gene not in genes:
                gene = 'intergenic'
            if gene not in mutations_dna:
                mutations_dna[gene] = []
                mutations_aa[gene] = []
            mutations_dna[gene].append(dna_name)
            if aa_name != 'nan' and gene != 'intergenic':
                mutations_aa[gene].append(aa_name)

        kmer_profile = list(set(kmer_profile_st))
        kmer_profile.sort()
        kmer_profile_st = ', '.join([str(x) for x in kmer_profile])
        kmer_profile_md5 = calc_md5(kmer_profile_st)
        genotype = ', '.join(list(data.keys()))
        mutations = []
        for gene in genes:
            if gene in mutations_dna:
                mutations.append(', '.join(list(set(mutations_dna[gene]))))
            else:
                mutations.append('')

        for gene in genes:
            if gene == 'intergenic':
                continue
            if gene in mutations_aa:
                mutations.append(', '.join(list(set(mutations_aa[gene]))))
            else:
                mutations.append('')

        num_targets_detected = len(kmer_summary[sample_id]['targets'])
        if num_targets_detected > 0:
            average_target_freq = sum(kmer_summary[sample_id]['counts']) / num_targets_detected
            target_freq_stdev = statistics.pstdev(kmer_summary[sample_id]['counts'])
        else:
            average_target_freq = 0
            target_freq_stdev = 0

        num_mixed_sites = len(kmer_summary[sample_id]['mixed'])

        if num_mixed_sites > max_mixed_sites:
            sample_type = 'multi'
        else:
            sample_type = 'single'

        primary_genotype = 'n/a'
        primary_genotype_transition_ratio = 'n/a'
        primary_genotype_slope = 'n/a'

        if sample_id in sample_complexity:
            if len(sample_complexity[sample_id]['primary_genotype']) > 0:
                primary_genotype = sample_complexity[sample_id]['primary_genotype']
                primary_genotype_transition_ratio = abs(sample_complexity[sample_id]['primary_genesis_ratio'])
                if len(sample_complexity[sample_id]['transision_slopes']) > 0:
                    primary_genotype_slope = abs(sample_complexity[sample_id]['transision_slopes'][0])

        if kmer_profile_md5 not in profile_names:
            phrase = ' '.join(generate_random_phrase())
            profile_names[kmer_profile_md5] = phrase
        else:
            phrase = profile_names[kmer_profile_md5]


        report.append("{}".format("\t".join([sample_id,
                                                genotype,
                                                sample_type,
                                                primary_genotype,
                                                str(primary_genotype_transition_ratio),
                                                str(primary_genotype_slope),
                                                str(num_targets_detected),
                                                str(average_target_freq),
                                                str(target_freq_stdev),
                                                str(num_mixed_sites),
                                                kmer_profile_st,
                                                kmer_profile_md5,
                                                phrase,
                                                quality['qc_sample_type'],
                                                quality['qc_max_missing_sites'],
                                                quality['qc_max_het_sites'],
                                                quality['qc_cov_stdev'],
                                                quality['qc_overall'],
                                                "\t".join(mutations)])))
    fh = open(outfile, 'w')
    fh.write("\n".join(report))


def perform_sample_complexity_analysis(sample_report, sample_kmer_data, min_cov_frac=0.05, step_size=0.05, min_cov=20):
    '''
    Using the sample kmer data and detected genotypes, it will determine the slope of the sample composition between each
    point where the number of samples differs
    :param sample_kmer_data: dict
    :param max_cov_frac: max
    :param min_cov_frac:
    :return:
    '''
    complexity = {}
    for sample_id in sample_report:
        complexity[sample_id] = {
            'genotypes': [],
            'x': [],
            'y': [],
            'slopes': [],
            'transision_points': [],
            'transision_slopes': [],
            'primary_genotype': '',
            'primary_genesis_ratio': -1
        }
        counts = sample_kmer_data[sample_id]['counts']
        ratios = sample_kmer_data[sample_id]['ratios']

        # determine the number of genotypes present at each step
        x = 1
        while x >= min_cov_frac:
            num_genotypes = 0
            genotypes = []
            count_valid = 0
            for genotype in sample_report[sample_id]:
                targets = sample_report[sample_id][genotype]['targets']
                count_valid = 0
                for target in targets:
                    if counts[target]['positive'] >= min_cov and ratios[target] >= x:
                        count_valid += 1
                if count_valid > 0:
                    num_genotypes += 1
                    genotypes.append(genotype)
            complexity[sample_id]['genotypes'].append(genotypes)
            complexity[sample_id]['x'].append(x)
            complexity[sample_id]['y'].append(num_genotypes)
            x -= step_size

        # calculate slope between each adjacent pair
        x = complexity[sample_id]['x']
        y = complexity[sample_id]['y']
        for i in range(0, len(x) - 1):
            k = i + 1
            complexity[sample_id]['slopes'].append((y[k] - y[i]) / (x[k] - x[i]))
            if y[i] == 1 and x[i] >= complexity[sample_id]['primary_genesis_ratio']:
                complexity[sample_id]['primary_genesis_ratio'] = x[i]
                complexity[sample_id]['primary_genotype'] = complexity[sample_id]['genotypes'][i][0]

        # determine the transition points
        slopes = complexity[sample_id]['slopes']

        for i in range(0, len(slopes) - 1):
            k = i + 1
            if slopes[i] != slopes[k]:
                x1 = x[i + 1]
                y1 = y[i + 1]
                x2 = x[k + 1]
                y2 = y[k + 1]
                complexity[sample_id]['transision_points'].append(x2)
                complexity[sample_id]['transision_slopes'].append((y2 - y1) / (x2 - x1))

    return complexity



def process_directory(data_dir, biohansel_fasta_file, outdir, logger, prefix, scheme_kmer_groups,
                      scheme_target_to_group_mapping,
                      scheme_kmer_target_info, min_cov, min_cov_frac, nthreads, batch_size=100):
    kmer_file = os.path.join(outdir, "{}.sample.bh.kmer.txt".format(prefix))
    summary_file = os.path.join(outdir, "{}.sample.bh.summary.txt".format(prefix))
    simple_file = os.path.join(outdir, "{}.sample.bh.simple.txt".format(prefix))
    logger.info("Identifying kmers which are found in input {}".format(data_dir))
    seq_files = get_sequence_files(data_dir)
    sample_kmer_results = {}
    for ftype in seq_files:

        num_files = len(seq_files[ftype])

        if num_files > 0 and num_files <= batch_size:
            (stdout, stderr) = bio_hansel.run_biohansel_directory(biohansel_fasta_file, data_dir, kmer_file,
                                                                  summary_file, simple_file, min_cov, min_cov_frac,
                                                                  nthreads)
            logger.info("Biohansel stdout: {}".format(stdout))
            logger.info("Biohansel sterr: {}".format(stderr))
            sample_kmer_biohansel_df = read_tsv(kmer_file)
            logger.info("Processing bioHansel results")
            sample_kmer_results.update(process_biohansel_kmer(scheme_kmer_groups, scheme_target_to_group_mapping,
                                                              scheme_kmer_target_info, sample_kmer_biohansel_df,
                                                              min_cov))

        elif num_files > batch_size:

            tracker = 0
            batch = []
            for i in range(0, num_files):
                tracker += 1
                batch.append(seq_files[ftype][i])
                if tracker >= batch_size or i == num_files - 1:
                    (stdout, stderr) = bio_hansel.run_biohansel_single(biohansel_fasta_file, batch, kmer_file,
                                                                       summary_file,
                                                                       simple_file, min_cov, min_cov_frac, nthreads)
                    tracker = 0
                    batch = []
                    sample_kmer_biohansel_df = read_tsv(kmer_file)
                    logger.info("Processing bioHansel results")
                    sample_kmer_results.update(
                        process_biohansel_kmer(scheme_kmer_groups, scheme_target_to_group_mapping,
                                               scheme_kmer_target_info, sample_kmer_biohansel_df, min_cov))
                    logger.info("Biohansel stdout: {}".format(stdout))
                    logger.info("Biohansel sterr: {}".format(stderr))
    return sample_kmer_results

def sample_qc(sample_kmer_data,total_targets,min_cov=20,min_cov_frac=0.05,max_mixed_sites=10, max_missing_sites=1400, sample_type='single',sample_cov_stdev_perc=0.75):
    kmer_summary = get_detected_target_summary(sample_kmer_data, min_cov=min_cov, min_cov_frac=min_cov_frac)
    for sample in sample_kmer_data:
        sample_kmer_data[sample]['quality'] = {
            'qc_sample_type':'Pass',
            'qc_max_missing_sites': 'Pass',
            'qc_max_het_sites': 'Pass',
            'qc_cov_stdev': 'Pass',
            'qc_messages': [],
            'qc_overall': 'Pass',
        }
        num_targets_detected = len(kmer_summary[sample]['targets'])
        missing_targets = num_targets_detected - total_targets
        if missing_targets > max_missing_sites:
            sample_kmer_data[sample]['quality']['qc_max_missing_sites'] = 'Fail'
            sample_kmer_data[sample]['quality']['qc_overall'] = 'Fail'

        mixed_sites = len(kmer_summary[sample]['mixed'])
        if mixed_sites > max_mixed_sites:
            sample_kmer_data[sample]['quality']['qc_max_het_sites'] = 'Fail'
            sample_kmer_data[sample]['quality']['qc_overall'] = 'Fail'
            if sample_type == 'single':
                sample_kmer_data[sample]['quality']['qc_sample_type'] = 'Fail'

        num_targets_detected = len(kmer_summary[sample]['targets'])
        if num_targets_detected > 0:
            average_target_freq = sum(kmer_summary[sample]['counts']) / num_targets_detected
            target_freq_stdev = statistics.pstdev(kmer_summary[sample]['counts'])
        else:
            average_target_freq = 0
            target_freq_stdev = 0

        if target_freq_stdev >= (average_target_freq*sample_cov_stdev_perc):
            sample_kmer_data[sample]['quality']['qc_cov_stdev'] = 'Fail'
            sample_kmer_data[sample]['quality']['qc_overall'] = 'Fail'


    return sample_kmer_data

def write_kmer_dict(kmer_data,out_file):
    fh = open(out_file,'w')
    fh.write("sample_id\ttarget\tseq\tis_positive\tfreq\n")
    for sample in kmer_data:
        for target in kmer_data[sample]:
            for seq in kmer_data[sample][target]:
                is_pos_kmer = str(kmer_data[sample][target][seq]['is_pos_kmer'])
                freq = str(kmer_data[sample][target][seq]['freq'])
                fh.write("{}\n".format("\t".join([sample,target,seq,is_pos_kmer,freq])))
    fh.close()


def identify_nearest_neighbors(sample_profiles,genotype_profile_data,max_dist=0.5,max_neighbors=10):
    profiles = dist_compatible_profiles(genotype_profile_data)
    data = {}
    for sample_id,profile in sample_profiles:
        profile = set(profile.split(','))
        n = {}
        for id in profiles:
            profile_2 = profiles[id]['profile']
            if profile == profile_2:
                j = 0
            if len(profile) > 0 or len(profile_2) >0:
                j = 1 - len(profile & profile_2) / len(profile | profile_2)
            else:
                j = 1
            if j <= max_dist:
                n[id] = j
        n = sorted(n.items(), key=lambda x: x[1], reverse=False)
        count = 0
        neighbors = {}
        for id in n:
            if count < max_neighbors:
                neighbors[id] = n[id]
            else:
                break
            count+=1
        data[sample_id] = neighbors
    return data




def bin_scheme_targets(scheme_df,max_length,window_size=500):
    bins = []
    mapping = {}
    for i in range(window_size,max_length+window_size,window_size):
        bins.append(i)

    for row in scheme_df.itertuples():
        target = str(row.index)
        position = int(row.kmer_start)

        for i in range(0,len(bins)):
            bin = bins[i]
            if position <= bin:
                mapping[target] = bin
                break

    return mapping

def write_genotype_report(sample_kmer_data,scheme_df,scheme_kmer_target_info,outfile,min_cov,min_cov_frac):
    genotype_targets = {}
    genotypes = []
    global_pos_targets = []
    for row in scheme_df.itertuples():
        target = row.key
        if isinstance(row.positive_seqs,float):
            positive_seqs = []
        else:
            positive_seqs = str(row.positive_seqs).split(',')
        if isinstance(row.partial_positive_seqs,float):
            partial_pos_seqs = []
        else:
            partial_pos_seqs = str(row.partial_positive_seqs).split(',')

        for g in positive_seqs:
            if g not in genotype_targets:
                genotype_targets[g] = {
                    'positive':[],
                    'partial':[],
                    'total':[]
                }
            global_pos_targets.append(str(target))
            genotype_targets[g]['positive'].append(str(target))
            genotype_targets[g]['total'].append(str(target))
            genotypes.append(g)

        for g in partial_pos_seqs:
            if g not in genotype_targets:
                genotype_targets[g] = {
                    'positive': [],
                    'partial': [],
                    'total': []
                }
            genotype_targets[g]['total'].append(str(target))
            genotype_targets[g]['partial'].append(str(target))
            genotypes.append(g)
    genotypes = list(set(genotypes))
    global_pos_targets = list(set(global_pos_targets))


    genotype_report = ["sample_id\tgenotype\tnum_scheme_targets_total\tnum_scheme_targets_positive\t"
                       "num_scheme_targets_partial\t"
                       "num_found_targets_positive\t"
                       "num_found_targets_partial\t"
                       "found_pos_mutations_dna\tfound_pos_mutations_aa\t"
                       "missing_pos_mutations_dna\tmissing_pos_mutations_aa\t"
                       "found_par_mutations_dna\tfound_par_mutations_aa\t"
                       "positive_exclusion_target\tnegative_exclusion_target"]
    for sample_id in sample_kmer_data:
        identified_positive_targets = []
        for target in global_pos_targets:
            cov = sample_kmer_data[sample_id]['counts'][target]['positive']
            ratio = sample_kmer_data[sample_id]['ratios'][target]
            if cov >= min_cov and ratio >= min_cov_frac:
                identified_positive_targets.append(target)
        print(identified_positive_targets)

        genotype_data = sample_kmer_data[sample_id]['genotypes']
        compatible_genotypes = list(set(genotype_data['include']) )
        for genotype in genotypes:
            if not genotype in genotype_targets:
                continue
            scheme_targets = genotype_targets[genotype]['total']

            found_targets = []
            scheme_pos_targets = genotype_targets[genotype]['positive']
            scheme_par_targets = genotype_targets[genotype]['partial']
            pos_exclusion = []
            neg_exclusion = []
            for target in scheme_targets:
                target = str(target)
                cov = sample_kmer_data[sample_id]['counts'][target]['positive'] + sample_kmer_data[sample_id]['counts'][target]['negative']
                ratio = sample_kmer_data[sample_id]['ratios'][target]
                if cov >= min_cov and ratio >= min_cov_frac:
                    found_targets.append(target)
                    if target not in scheme_pos_targets:
                        neg_exclusion.append(target)
                elif cov >= min_cov:
                    if target in scheme_pos_targets:
                        pos_exclusion.append(target)

            found_pos_targets = list(set(scheme_pos_targets) & set(found_targets))
            found_par_targets = list(set(scheme_par_targets) & set(found_targets))
            missing_pos_targets = list(set(scheme_pos_targets) - set(found_targets))
            found_pos_mut_dna = []
            found_pos_mut_aa = []
            missing_pos_mut_dna = []
            missing_pos_mut_aa = []
            found_par_mut_dna = []
            found_par_mut_aa = []
            for target in found_pos_targets:
                found_pos_mut_dna.append(str(scheme_kmer_target_info[target]['dna_name']))
                aa = str(scheme_kmer_target_info[target]['aa_name'])
                if aa == 0 or aa == '0' or aa == '':
                    aa = "N/A_{}".format(str(scheme_kmer_target_info[target]['dna_name']))
                found_pos_mut_aa.append(str(aa))
            for target in missing_pos_targets:
                missing_pos_mut_dna.append(str(scheme_kmer_target_info[target]['dna_name']))
                aa = str(scheme_kmer_target_info[target]['aa_name'])
                if aa == 0 or aa == '0' or aa == '':
                    aa = "N/A_{}".format(str(scheme_kmer_target_info[target]['dna_name']))
                missing_pos_mut_aa.append(str(aa))
            for target in found_par_targets:
                found_par_mut_dna.append(str(scheme_kmer_target_info[target]['dna_name']))
                aa = str(scheme_kmer_target_info[target]['aa_name'])
                if aa == 0 or aa == '0' or aa == '':
                    aa = "N/A_{}".format(str(scheme_kmer_target_info[target]['dna_name']))
                found_par_mut_aa.append(str(aa))
            pos_exclusion_mut = []
            for target in pos_exclusion:
                pos_exclusion_mut.append(str(scheme_kmer_target_info[target]['dna_name']))
            neg_exclusion_mut = []
            for target in neg_exclusion:
                neg_exclusion_mut.append(str(scheme_kmer_target_info[target]['dna_name']))

            row = [sample_id,genotype,str(len(scheme_targets)),
                   str(len(scheme_pos_targets)),
                   str(len(scheme_par_targets)),
                   str(len(found_pos_targets)),
                   str(len(found_par_targets)),
                   ", ".join(found_pos_mut_dna),
                   ", ".join(found_pos_mut_aa),
                   ", ".join(missing_pos_mut_dna),
                   ", ".join(missing_pos_mut_aa),
                   ", ".join(found_par_mut_dna),
                   ", ".join(found_pos_mut_aa),
                   ", ".join(pos_exclusion_mut),
                   ", ".join(neg_exclusion_mut),
                   ]

            genotype_report.append("\t".join(row))
    #dna_name = str(scheme_kmer_target_info[target]['dna_name'])

    fh = open(outfile,'w')
    fh.write("\n".join(genotype_report))
    fh.close




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
    scheme = cmd_args.scheme
    mode = cmd_args.mode
    type = cmd_args.type
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
    report_run_info_log.write("Genotype Profiles\t{}\n".format(genotype_profiles))
    report_run_info_log.write("Positive Control\t{}\n".format(positive_control))
    report_run_info_log.write("No Template Control\t{}\n".format(no_template_control))
    report_run_info_log.write("Primers\t{}\n".format(primers))

    logger.info("Reading kmer scheme from {}".format(scheme))
    if scheme in TYPING_SCHEMES:
        scheme = TYPING_SCHEMES[scheme]
    report_run_info_log.write("Scheme\t{}\n".format(scheme))
    scheme_df = read_tsv(scheme)
    genes = scheme_df['gene'].dropna().unique().tolist()
    genes.append('intergenic')
    logger.info("Found {} kmers in {}".format(len(scheme_df), scheme))
    report_run_info_log.write("Number of scheme kmer targets\t{}\n".format(len(scheme_df)))
    biohansel_fasta_file = os.path.join(outdir, "{}.biohansel.fasta".format(prefix))
    logger.info("Writing biohansel compatible kmer scheme to {}".format(biohansel_fasta_file))
    scheme_bin_mapping = bin_scheme_targets(scheme_df,30000,window_size=500)
    scheme_to_biohansel_fasta(scheme_df, biohansel_fasta_file)
    scheme_kmer_target_info = init_kmer_targets(scheme_df)
    scheme_kmer_target_keys = list(scheme_kmer_target_info.keys())
    scheme_kmer_target_keys.sort()
    scheme_kmer_target_keys = [str(i) for i in scheme_kmer_target_keys]
    scheme_kmer_groups = get_kmer_groups(scheme_kmer_target_info)
    scheme_target_to_group_mapping = get_kmer_group_mapping(scheme_kmer_target_info)

    # Identify kmers which are present in the primer sets for masking purposes later
    primer_kmers = {}
    if primers is not None:
        if primers in PRIMER_SCHEMES:
            primers = PRIMER_SCHEMES[primers]
        logger.info("Writing primers to fasta file".format(primers))
        primer_df = read_tsv(primers)
        primer_file = os.path.join(outdir, "primers.fasta")
        fh = open(primer_file, 'w')
        for row in primer_df.itertuples():
            fh.write(">{}\n{}\n".format(row.name, row.seq))
        fh.close()

        kmer_file = os.path.join(outdir, "primers.bh.kmer.txt")
        summary_file = os.path.join(outdir, "primers.bh.summary.txt")
        simple_file = os.path.join(outdir, "primers.bh.simple.txt")
        logger.info("Identifying kmers which are found in selected primer set {}".format(primers))
        (stdout, stderr) = bio_hansel.run_biohansel_single(biohansel_fasta_file, [primer_file], kmer_file, summary_file,
                                                           simple_file, min_cov, min_cov_frac, nthreads)
        logger.info("Biohansel stdout: {}".format(stdout))
        logger.info("Biohansel sterr: {}".format(stderr))
        primer_df = read_tsv(kmer_file)
        primer_data = process_biohansel_kmer(scheme_kmer_groups, scheme_target_to_group_mapping,
                                             scheme_kmer_target_info, primer_df, min_cov)

        write_kmer_dict(primer_data, report_primers)
        del (primer_df)
        primer_kmers = {}
        if len(primer_data) > 0:
            for target in primer_data['primers']:
                for seq in primer_data['primers'][target]:
                    if primer_data['primers'][target][seq]['is_pos_kmer']:
                        kmername = target
                    else:
                        kmername = "negative{}".format(target)
                    primer_kmers[kmername] = min_cov

            report_run_info_log.write("Number of kmers found in primers\t{}\n".format(len(primer_kmers)))
        logger.info("Found {} kmers overlapping primers".format(len(primer_kmers)))

    # Identify kmers which are present in the no template control for masking purposes later
    no_template_kmers = {}
    no_template_results = {}
    if no_template_control is not None:
        kmer_file = os.path.join(outdir, "{}.no_template.bh.kmer.txt".format(prefix))
        summary_file = os.path.join(outdir, "{}.no_template.bh.summary.txt".format(prefix))
        simple_file = os.path.join(outdir, "{}.no_template.bh.simple.txt".format(prefix))
        logger.info("Identifying kmers which are found in no template control {}".format(no_template_control))
        (stdout, stderr) = bio_hansel.run_biohansel_single(biohansel_fasta_file, [no_template_control], kmer_file,
                                                           summary_file, simple_file, min_cov, nthreads)
        logger.info("Biohansel stdout: {}".format(stdout))
        logger.info("Biohansel sterr: {}".format(stderr))
        no_template_df = read_tsv(kmer_file)
        no_template_results = process_biohansel_kmer(scheme_kmer_groups, scheme_target_to_group_mapping,
                                                     scheme_kmer_target_info, no_template_df, min_cov)

        if len(no_template_results) > 0:
            sample = list(no_template_results.keys())[0]
            write_kmer_dict(no_template_results, report_negative)
            for kmer in no_template_results[sample]:
                for seq in no_template_results[sample][kmer]:
                    freq = no_template_results[sample][kmer][seq]['freq']
                    if no_template_results[sample][kmer][seq]['is_pos_kmer']:
                        kmername = kmer
                    else:
                        kmername = "negative{}".format(kmer)
                    if kmername not in no_template_kmers:
                        no_template_kmers[kmername] = 0
                    no_template_kmers[kmername] += freq
            report_run_info_log.write("No template kmer targets\t{}\n".format(len(no_template_kmers)))
            no_template_results = calc_kmer_ratio(no_template_results, scheme_kmer_target_keys, min_cov)
        logger.info("Found {} kmers found in no template".format(len(no_template_kmers)))


    report_run_info_log.write("No Template Control kmers\t{}\n".format(len(no_template_kmers)))
    # Identify kmers identified in positive control
    positive_control_kmers = {}
    positive_control_results = {}
    if positive_control is not None:
        kmer_file = os.path.join(outdir, "positive_control.bh.kmer.txt")
        summary_file = os.path.join(outdir, "positive_control.bh.summary.txt")
        simple_file = os.path.join(outdir, "positive_control.bh.simple.txt")
        logger.info("Identifying kmers which are found in positive control {}".format(positive_control))
        (stdout, stderr) = bio_hansel.run_biohansel_single(biohansel_fasta_file, [positive_control], kmer_file,
                                                           summary_file, simple_file, min_cov, min_cov_frac, nthreads)
        positive_control_df = read_tsv(kmer_file)
        positive_control_results = process_biohansel_kmer(scheme_kmer_groups, scheme_target_to_group_mapping,
                                                          scheme_kmer_target_info, positive_control_df, min_cov)
        write_kmer_dict(positive_control_results, report_positive)

        if len(positive_control_results) > 0:
            sample = list(positive_control_results.keys())[0]
            num_kmers_found = len(positive_control_results[sample])
            for kmer in positive_control_results[sample]:
                for seq in positive_control_results[sample][kmer]:
                    freq = positive_control_results[sample][kmer][seq]['freq']
                    if positive_control_results[sample][kmer][seq]['is_pos_kmer']:
                        kmername = kmer
                    else:
                        kmername = "negative{}".format(kmer)
                    if kmername not in positive_control_kmers:
                        positive_control_kmers[kmername] = 0
                    positive_control_kmers[kmername] += freq
            positive_control_results = calc_kmer_ratio(positive_control_results, scheme_kmer_target_keys, min_cov)
        logger.info("Found {} kmers in positive control".format(len(positive_control_results[sample])))
        report_run_info_log.write("Positive countrol kmer targets\t{}\n".format(num_kmers_found))
        logger.warning("Missing {} kmers in positive control".format(len(scheme_df) - num_kmers_found))
        report_run_info_log.write(
            "Missing Positive countrol kmer targets\t{}\n".format(len(scheme_df) - len(positive_control_kmers)))
        if (len(scheme_df) - len(positive_control_kmers)) > max_missing_sites:
            check = 'Fail'
            logger.warning("Positive control has lower than required number of k-mer targets")
        else:
            check = 'Pass'
        report_run_info_log.write(
            "Positive countrol QC Check\t{}\n".format(check))

        logger.info("Biohansel stdout: {}".format(stdout))
        logger.info("Biohansel sterr: {}".format(stderr))

    kmer_file = os.path.join(outdir, "{}.sample.bh.kmer.txt".format(prefix))
    summary_file = os.path.join(outdir, "{}.sample.bh.summary.txt".format(prefix))
    simple_file = os.path.join(outdir, "{}.sample.bh.simple.txt".format(prefix))

    # Run kmer detection on the samples
    if data_dir is not None:
        input_data = data_dir
        logger.info("Identifying kmers which are found in input directory {}".format(data_dir))
        logger.info("Processing bioHansel results")
        sample_kmer_results = process_directory(data_dir, biohansel_fasta_file, outdir, logger, prefix,
                                                scheme_kmer_groups,
                                                scheme_target_to_group_mapping,
                                                scheme_kmer_target_info, min_cov, min_cov_frac, nthreads,
                                                batch_size=100)

    elif SE is not None:
        input_data = SE

        logger.info("Identifying kmers which are found in input {}".format(SE))
        (stdout, stderr) = bio_hansel.run_biohansel_single(biohansel_fasta_file, [SE], kmer_file, summary_file,
                                                           simple_file, min_cov, min_cov_frac, nthreads)
        sample_kmer_biohansel_df = read_tsv(kmer_file)
        logger.info("Processing bioHansel results")
        sample_kmer_results = process_biohansel_kmer(scheme_kmer_groups, scheme_target_to_group_mapping,
                                                     scheme_kmer_target_info, sample_kmer_biohansel_df, min_cov)
    else:
        input_data = "{},{}".format(R1, R2)
        logger.info("Identifying kmers which are found in input {} {}".format(R1, R2))
        (stdout, stderr) = bio_hansel.run_biohansel_pair(biohansel_fasta_file, R1, R2, kmer_file, summary_file,
                                                         simple_file, min_cov, min_cov_frac, nthreads)
        sample_kmer_biohansel_df = read_tsv(kmer_file)
        logger.info("Processing bioHansel results")
        sample_kmer_results = process_biohansel_kmer(scheme_kmer_groups, scheme_target_to_group_mapping,
                                                     scheme_kmer_target_info, sample_kmer_biohansel_df, min_cov)

#    logger.info("Biohansel stdout: {}".format(stdout))
 #   logger.info("Biohansel sterr: {}".format(stderr))
    logger.info("Reading bioHansel results")
    report_run_info_log.write("Input Data\t{}\n".format(input_data))
    logger.info("Determining kmer target ratios")
    sample_kmer_data = calc_kmer_ratio(sample_kmer_results, scheme_kmer_target_keys, min_cov)
    logger.info("Identifying mixed kmer targets")
    sample_kmer_data = calc_mixed_sites(sample_kmer_data, min_cov_frac)
    report_run_info_log.write("Number of samples\t{}\n".format(len(sample_kmer_data)))
    report_run_info_log.write("Positive control kmers\t{}\n".format(len(positive_control_kmers)))

    # subtract contamination kmers in the no template control
    if len(no_template_kmers) > 0:
        logger.info("Subtracting no template control kmers from input data")
        sample_kmer_data = subtract_contamination_kmers(sample_kmer_data, no_template_kmers, min_cov_frac)

    # subtract potential contamination kmers from the primers
    if len(primer_kmers) > 0:
        ##To do, be able to distinguish where in the read the kmer is to address the primer issue
        logger.warn("Warning found kmers overlapping primers, recommend preprocessing reads to remove primer sequences")
        sample_kmer_data = subtract_contamination_kmers(sample_kmer_data, primer_kmers, min_cov_frac)

    sample_kmer_data = calc_mixed_sites(sample_kmer_data, min_cov_frac)
    # Get information about each of the targets involved in mixed events
    target_info = get_mixed_kmer_targets(sample_kmer_data, min_cov_frac)
    target_report = pd.DataFrame().from_dict(target_info, orient='index')
    target_report.reset_index(level=0, inplace=True)
    target_report = target_report.rename(columns={'index': 'kmer_target'})
    target_report.to_csv(report_run_kmer_composition, sep="\t", index=False)

    # Get list of strains which are compatible with the kmer information
    sample_kmer_data = identify_compatible_types(scheme_df, sample_kmer_data, min_cov_frac,
                                                 detection_limit=detection_limit)

    sample_kmer_data = calc_type_coverage(sample_kmer_data, scheme_df, min_cov_frac=min_cov_frac, min_cov=min_cov)
    sample_kmer_data = type_occamization(sample_kmer_data, scheme_df, min_cov_frac=min_cov_frac, min_cov=min_cov)
    sample_kmer_data = calc_type_coverage(sample_kmer_data, scheme_df, min_cov_frac=min_cov_frac, min_cov=min_cov,
                                          recalc=True)
    profile_st = {}
    logger.info("Calculating scipy compatible profile features")
    kmer_content_profile_df = generate_target_presence_table(sample_kmer_data, min_ratio=min_cov_frac, max_ratio=1-min_cov_frac)
    kmer_content_profile_df['target'] = kmer_content_profile_df.index
    first_column = kmer_content_profile_df.pop('target')
    kmer_content_profile_df.insert(0, 'target', first_column)
    logger.info("Writing scipy compatible profile to: {}".format(report_scipy_profile))
    kmer_content_profile_df.to_csv(report_scipy_profile,header=True,sep='\t',index=False)
    if len(sample_kmer_data) > 1:
        plot_mds(nan_compatible_kmer_pairwise_distmatrix(kmer_content_profile_df), list(kmer_content_profile_df.columns),report_run_kmer_mds)

    # get positive kmer signature for each sample
    kmer_summary = get_detected_target_summary(sample_kmer_data, min_cov, min_cov_frac)

    sample_kmer_data = sample_qc(sample_kmer_data, len(scheme_df), min_cov=min_cov, min_cov_frac=min_cov_frac, max_mixed_sites=max_mixed_sites,
              max_missing_sites=max_missing_sites, sample_type=type, sample_cov_stdev_perc=0.75)

    for sample_id in kmer_summary:
        profile_st[sample_id] = kmer_summary[sample_id]['positive']
    genotype_profiles_data = {}
    neighbours = {}
    if genotype_profiles is not None:
        genotype_profiles_data = read_genotype_profiles(genotype_profiles)
        neighbours = identify_nearest_neighbors(profile_st,genotype_profiles_data,max_dist=genotype_dist_cutoff,max_neighbors=10)

    sample_report = {}
    fail_sample_count = 0
    for sample in sample_kmer_data:
        if sample_kmer_data[sample]['quality']['qc_overall'] == 'Fail':
            fail_sample_count+=1
        if not sample in sample_report:
            sample_report[sample] = {}
        total_ratio = 0
        unassigned_positive_kmers = kmer_summary[sample]['positive']

        for genotype in sample_kmer_data[sample]['genotypes']['candidate_data']:
            data = sample_kmer_data[sample]['genotypes']['candidate_data'][genotype]

            # Filter out strains with low confidence
            if (data['average_ratio'] > 0 and data['average_ratio'] < min_cov_frac) or data[
                'average_target_freq'] < min_cov:
                continue

            data['targets'].sort()
            unassigned_positive_kmers = list(set(unassigned_positive_kmers) - set(data['targets']))
            total_ratio += data['average_ratio']
            sample_report[sample][genotype] = {
                'num_targets': data['num_targets'],
                'targets': data['targets'],
                'average_target_freq': data['average_target_freq'],
                'target_freq_stdev': data['target_freq_stdev'],
                'average_ratio': data['average_ratio'],
                'ratio_stdev': data['ratio_stdev']
            }

        # Assign kmers not assigned to a genotype to unknown genotype
        unknown_fraq = 1 - total_ratio
        if len(unassigned_positive_kmers) > 0 and unknown_fraq > min_cov_frac:
            counts = []
            ratios = []
            for kmer_id in unassigned_positive_kmers:
                if kmer_id in sample_kmer_data[sample]['counts']:
                    counts.append(sample_kmer_data[sample]['counts'][kmer_id]['positive'])
                    ratios.append(sample_kmer_data[sample]['ratios'][kmer_id])

            ave_count = sum(counts) / len(counts)
            std_count = statistics.pstdev(counts)
            ave_ratio = sum(ratios) / len(ratios)
            std_ratio = statistics.pstdev(counts)

            sample_report[sample]['unknown'] = {
                'num_targets': len(unassigned_positive_kmers),
                'targets': unassigned_positive_kmers,
                'average_target_freq': ave_count,
                'target_freq_stdev': std_count,
                'average_ratio': ave_ratio,
                'ratio_stdev': std_ratio
            }

    sample_complexity = perform_sample_complexity_analysis(sample_report, sample_kmer_data, min_cov_frac=min_cov_frac,
                                                           min_cov=min_cov, step_size=0.05)
    write_sample_detailed_report(sample_report, report_sample_composition_detailed, scheme_kmer_target_info)
    write_sample_summary_report(sample_report, report_sample_composition_summary, sample_kmer_data,
                                scheme_kmer_target_info, sample_complexity, genes, max_mixed_sites=max_mixed_sites,
                                min_cov_frac=min_cov_frac, min_cov=min_cov,sample_type=type)
    write_genotype_report(sample_kmer_data, scheme_df, scheme_kmer_target_info, report_genotype_targets, min_cov, min_cov_frac)
    # create a plot of sample similarity for a multi-sample run
    if len(profile_st) > 1:
        labels = []
        for sample in profile_st:
            genotype = ', '.join(list(sample_report[sample].keys()))
            labels.append("{} | {}".format(sample,genotype))
        d = dendrogram_visualization()
        d.build_tree_from_dist_matrix(labels, profile_pairwise_distmatrix(profile_st),
                                      report_run_kmer_dendrogram)
    report_run_info_log.write("Samples failing QC\t{}\n".format(fail_sample_count))
    report_run_info_log.write("End Time\t{}\n".format(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
    coverage_summary_fig = generate_sample_coverage_plot(positive_control_results, no_template_results, sample_kmer_data, scheme_bin_mapping)
    coverage_summary_fig.write_image(report_run_kmer_coverage)
