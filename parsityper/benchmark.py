import sys, os
from argparse import (ArgumentParser, FileType)
from itertools import permutations
from parsityper.kmerSearch import init_automaton_dict,find_in_fasta_dict
from parsityper.helpers import read_tsv, read_fasta, init_kmer_targets, get_kmer_groups, \
    get_kmer_group_mapping, process_biohansel_kmer
from parsityper.typer import calc_kmer_ratio, calc_mixed_sites, identify_compatible_types, \
calc_type_coverage, type_occamization, get_detected_target_summary, get_genotype_target_associations
import statistics
import time
from multiprocessing import Pool
from multiprocessing import set_start_method
set_start_method("spawn")

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Kmer scheme generator')
    parser.add_argument('--input_seqs', type=str, required=True,
                        help='single fasta format with all seqs for testing')
    parser.add_argument('--outdir', type=str, required=True,
                        help='output directory')
    parser.add_argument('--input_meta', type=str, required=False,
                        help='TSV file of sample_id,genotype')
    parser.add_argument('--scheme', type=str, required=True,
                        help='TSV formated kmer scheme', default='')
    parser.add_argument('--n_threads', type=int, required=False,
                        help='Number of threads to use',default=1)

    return parser.parse_args()

def generate_comparsions(samples,depth=2):
    return list(permutations(samples, depth))

def call_genotype(sample_kmer_data,kmer_summary,combo,min_cov,min_cov_frac):
    sample_report = {}
    for sample in sample_kmer_data:
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
    overlap = list(set(combo) & set(list(sample_report[sample].keys())))
    return ("{}\t{}\t{}\n".format(sample, ",".join(list(sample_report[sample].keys())), len(overlap)))

def process_sample(sample_id,combo,sample_kmer_results_main,scheme_kmer_target_keys,scheme_df,min_cov,min_cov_frac):
    sample_kmer_results = {}
    sample_kmer_results[sample_id] = {}
    for cid in combo:
        for t in sample_kmer_results_main[cid]:
            if t not in sample_kmer_results[sample_id]:
                sample_kmer_results[sample_id][t] = {}
            for seq in sample_kmer_results_main[cid][t]:
                if not seq in sample_kmer_results[sample_id][t]:
                    sample_kmer_results[sample_id][t][seq] = {'freq': 0,
                                                              'is_pos_kmer': sample_kmer_results_main[cid][t][seq][
                                                                  'is_pos_kmer']}
                sample_kmer_results[sample_id][t][seq]['freq'] += sample_kmer_results_main[cid][t][seq]['freq']

    sample_kmer_data = calc_kmer_ratio(sample_kmer_results, scheme_kmer_target_keys, 20)
    sample_kmer_data = calc_mixed_sites(sample_kmer_data, min_cov_frac)

    # Get list of strains which are compatible with the kmer information
    sample_kmer_data = identify_compatible_types(scheme_df, sample_kmer_data, min_cov_frac,
                                                 detection_limit=10)

    sample_kmer_data = calc_type_coverage(sample_kmer_data, scheme_df, min_cov_frac=min_cov_frac, min_cov=20)
    sample_kmer_data = type_occamization(sample_kmer_data, scheme_df, min_cov_frac=min_cov_frac, min_cov=20)
    sample_kmer_data = calc_type_coverage(sample_kmer_data, scheme_df, min_cov_frac=min_cov_frac, min_cov=20,
                                          recalc=True)
    kmer_summary = get_detected_target_summary(sample_kmer_data, min_cov, min_cov_frac)
    return call_genotype(sample_kmer_data, kmer_summary, combo, min_cov, min_cov_frac)

def run():
    cmd_args = parse_args()
    scheme = cmd_args.scheme
    input_seqs = cmd_args.input_seqs
    outdir = cmd_args.outdir
    scheme_df = read_tsv(scheme)
    n_threads = cmd_args.n_threads

    min_cov_frac = 0.4
    min_cov = 20
    # initialize analysis directory
    if not os.path.isdir(outdir):
        os.mkdir(outdir, 0o755)

    out_file = os.path.join(outdir,'pariwise.txt')
    out = open(out_file,'w')
    scheme_dict = {}
    for row in scheme_df.itertuples():
        id = row.key
        pos_kmer = row.positive
        neg_kmer = row.negative
        scheme_dict["{}-{}".format(id,id)] = pos_kmer
        scheme_dict["negative{}-{}".format(id,id)] = neg_kmer

    scheme_kmer_target_info = init_kmer_targets(scheme_df)
    scheme_kmer_target_keys = list(scheme_kmer_target_info.keys())
    scheme_kmer_target_keys.sort()
    scheme_kmer_target_keys = [str(i) for i in scheme_kmer_target_keys]
    scheme_kmer_groups = get_kmer_groups(scheme_kmer_target_info)
    scheme_target_to_group_mapping = get_kmer_group_mapping(scheme_kmer_target_info)

    A = init_automaton_dict(scheme_dict)
    seqs = read_fasta(input_seqs)
    kmer_df = find_in_fasta_dict(A,seqs)
    refpositions = [x for x, y in kmer_df.kmername.str.split('-')]
    kmer_df['refposition'] = [int(x.replace('negative', '')) for x in refpositions]
    kmer_df['is_pos_kmer'] = ~kmer_df.kmername.str.contains('negative')
    kmer_df['sample'] = kmer_df['contig_id']
    #print(kmer_df.columns)
    sample_kmer_results_main = process_biohansel_kmer(scheme_kmer_groups, scheme_target_to_group_mapping,
                                                 scheme_kmer_target_info, kmer_df, 20)

    #print(sample_kmer_results_main)
    combinations = generate_comparsions(list(sample_kmer_results_main.keys()),depth=2)

    sample_ids = []
    sample_report = {}
    pool = Pool(processes=n_threads)
    res = []
    results = []
    for i in range(0,len(combinations)):
        stime = time.time()
        combo = list(combinations[i])
        combo.sort()
        sample_id = ':'.join(combo)

        if sample_id in sample_ids:
            continue
        else:
            sample_ids.append(sample_id)
        res.append(pool.apply_async(process_sample,(sample_id,combo,sample_kmer_results_main,scheme_kmer_target_keys,scheme_df,min_cov,min_cov_frac)))
        print(sample_id)
        sys.stdin.flush()
        if len(res) == 100:
            pool.close()
            pool.join()
            results.extend([x.get() for x in res])
            pool = Pool(processes=n_threads)
            res = []
            print("{}\t{}".format(time.time()-stime,len(results)))
            stime = time.time()

    for result in results:
        print("{}".format(result))
