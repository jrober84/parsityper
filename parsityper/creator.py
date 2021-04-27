#!/usr/bin/python
from argparse import (ArgumentParser, FileType)
import logging, os, sys, re, collections, operator, math, shutil, datetime
from collections import Counter
import pandas as pd
from parsityper.helpers import init_console_logger, read_tsv, parse_reference_sequence, read_fasta, calc_consensus, generate_consensus_seq,\
find_snp_positions, find_snp_kmers, find_indel_kmers, find_internal_gaps, scheme_to_biohansel_fasta, init_kmer_targets, count_kmers, count_ambig, \
get_kmer_complexity, calc_md5
from parsityper.bio_hansel import bio_hansel
from parsityper.helpers import  profile_pairwise_distmatrix
from parsityper.visualizations import dendrogram_visualization

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Kmer scheme generator')
    parser.add_argument('--input_msa', type=str, required=True,
                        help='MSA in fasta format')
    parser.add_argument('--input_meta', type=str, required=True,
                        help='TSV file of sample_id,genotype')
    parser.add_argument('--ref_id', type=str, required=True,
                        help='sample_id for reference sequence to use in MSA')
    parser.add_argument('--ref_gbk', type=str, required=True,
                        help='GenBank file for reference sequences')
    parser.add_argument('--outdir', type=str, required=True,
                        help='output directory')
    parser.add_argument('--prefix', type=str, required=False,
                        help='output file prefix',default='parsityper')
    parser.add_argument('--min_len', type=int, required=False,
                        help='Absolute minimum length of acceptable k-mer',default=17)
    parser.add_argument('--max_len', type=int, required=False,
                        help='Absolute minimum length of acceptable k-mer',default=21)
    parser.add_argument('--max_ambig', type=int, required=False,
                        help='Absolute maximum of degenerate bases allowed in a k-mer',default=4)
    parser.add_argument('--max_states', type=int, required=False,
                        help='Absolute maximum of states allowed per kmer',default=256)
    parser.add_argument('--min_complexity', type=int, required=False,
                        help='Absolute maximum of degenerate bases allowed in a k-mer',default=0.2)

    return parser.parse_args()

def print_scheme(scheme,file):
    """
    Takes the kmer dictionary and writes it to a file
    :param scheme: dict of kmer info
    :param file: file path
    :return:
    """
    df = pd.DataFrame().from_dict(scheme,orient='index')
    df.reset_index(inplace=True)
    df['key'] = df.index
    df.to_csv(file,sep="\t",index=False)

def build_kmer_profiles(scheme):
    '''
    :param scheme:
    :type scheme:
    :return:
    :rtype:
    '''
    profiles = {}
    for kmer_id in scheme:
        sequences = scheme[kmer_id]['positive_seqs']
        for sample_id in sequences:
            if not sample_id in profiles:
                profiles[sample_id] = []
            profiles[sample_id].append(kmer_id)
    return profiles

def get_genotype_mapping(metadata_df):
    mapping = {}
    for row in metadata_df.itertuples():
        sample_id = row.sample_id
        genotype = row.genotype
        mapping[sample_id] = genotype
    return mapping



def qa_scheme(scheme,missing_targets, multiplicity_targets,min_len,max_len,max_ambig,min_complexity):
    '''

    :param scheme:
    :type scheme:
    :param missing_targets:
    :type missing_targets:
    :param multiplicity_targets:
    :type multiplicity_targets:
    :param min_len:
    :type min_len:
    :param max_len:
    :type max_len:
    :param max_ambig:
    :type max_ambig:
    :param min_complexity:
    :type min_complexity:
    :return:
    :rtype:
    '''
    kmer_start = []
    for kmer_id in scheme:
        pos_kmer = scheme[kmer_id]['positive'].replace('-','')
        neg_kmer = scheme[kmer_id]['negative'].replace('-','')
        pos_kmer_len = len(pos_kmer)
        scheme[kmer_id]['pos_kmer_len'] = pos_kmer_len
        scheme[kmer_id]['pos_kmer_len'] = pos_kmer_len >= min_len and pos_kmer_len <= max_len
        neg_kmer_len = len(neg_kmer)
        scheme[kmer_id]['neg_kmer_len'] = neg_kmer_len
        scheme[kmer_id]['neg_kmer_len'] = neg_kmer_len >= min_len and neg_kmer_len <= max_len
        bases = ['A','T','C','G']
        pos_kmer = list(scheme[kmer_id]['positive'])
        neg_kmer = list(scheme[kmer_id]['negative'])

        #remove degenerate bases which are common and at the begining or end
        for i in range(0,pos_kmer_len):
            pbase = pos_kmer[i]
            nbase = neg_kmer[i]

            if pbase != nbase or \
                   pbase in bases or \
                    nbase in bases:
                break
            if len(''.join(pos_kmer).replace('-','')) <= min_len or len(''.join(neg_kmer).replace('-','')) <= min_len :
                break

            scheme[kmer_id]['kmer_start']-=1
            pos_kmer[i] = '-'
            neg_kmer[i] = '-'

        for i in range(pos_kmer_len-1,0):
            pbase = pos_kmer[i]
            nbase = neg_kmer[i]

            if pbase != nbase or \
                   pbase in bases or \
                    nbase in bases:
                break
            if len(''.join(pos_kmer).replace('-','')) <= min_len or len(''.join(neg_kmer).replace('-','')) <= min_len :
                break
            scheme[kmer_id]['kmer_end'] -= 1
            pos_kmer[i] = '-'
            neg_kmer[i] = '-'

        if not scheme[kmer_id]['kmer_start'] in kmer_start:
            kmer_start.append(scheme[kmer_id]['kmer_start'])
        else :
            i = 0
            while scheme[kmer_id]['kmer_start'] in kmer_start:
                if pos_kmer_len < min_len or pos_kmer[i] != neg_kmer[i]:
                    break
                scheme[kmer_id]['kmer_start'] -= 1
                pos_kmer[i] = '-'
                neg_kmer[i] = '-'
                pos_kmer_len = len(''.join(pos_kmer).replace('-', ''))
                i+=1
            kmer_start.append(scheme[kmer_id]['kmer_start'])


        pos_kmer = ''.join(pos_kmer).replace('-','')
        neg_kmer = ''.join(neg_kmer).replace('-','')
        scheme[kmer_id]['positive'] = pos_kmer
        scheme[kmer_id]['negative'] = neg_kmer
        base_counts = count_kmers(pos_kmer, 1)

        num_bases_found = 0
        sum_gc = 0
        for b in bases:
            if b in base_counts:
                num_bases_found+=1
                if b == 'G' or b == 'C':
                    sum_gc+= base_counts[b]
        perc_gc = sum_gc / sum(list(base_counts.values()))
        scheme[kmer_id]['gc'] = perc_gc

        pos_ambig_count = count_ambig(pos_kmer)
        neg_ambig_count = count_ambig(neg_kmer)
        pos_complexity = get_kmer_complexity(pos_kmer)
        neg_complexity = get_kmer_complexity(neg_kmer)
        scheme[kmer_id]['pos_ambig_count'] = pos_ambig_count
        scheme[kmer_id]['pos_ambig_count_ok'] = pos_ambig_count <= max_ambig
        scheme[kmer_id]['neg_ambig_count'] = neg_ambig_count
        scheme[kmer_id]['neg_ambig_count_ok'] = neg_ambig_count <= max_ambig
        scheme[kmer_id]['pos_complexity'] = pos_complexity['average']
        scheme[kmer_id]['pos_complexity_ok'] = pos_complexity['average'] <= min_complexity
        scheme[kmer_id]['neg_complexity'] = neg_complexity['average']
        scheme[kmer_id]['neg_complexity_ok'] = neg_complexity['average'] <= min_complexity
        scheme[kmer_id]['is_unique'] = True
        scheme[kmer_id]['is_detectable'] = True
        scheme[kmer_id]['is_valid'] = True

        if not scheme[kmer_id]['pos_kmer_len'] > min_len or scheme[kmer_id]['neg_kmer_len'] > min_len:
            scheme[kmer_id]['is_valid'] = False

        if not scheme[kmer_id]['pos_complexity_ok']:
            scheme[kmer_id]['is_valid'] = False
        if not scheme[kmer_id]['neg_complexity_ok']:
            scheme[kmer_id]['is_valid'] = False
        if not scheme[kmer_id]['pos_ambig_count_ok']:
            scheme[kmer_id]['is_valid'] = False
        if not scheme[kmer_id]['neg_ambig_count_ok']:
            scheme[kmer_id]['is_valid'] = False

        if kmer_id in multiplicity_targets:
            scheme[kmer_id]['is_unique'] = False
            scheme[kmer_id]['is_valid'] = False
        if kmer_id in missing_targets:
            scheme[kmer_id]['is_detectable'] = False
            scheme[kmer_id]['is_valid'] = False

    return scheme

def identify_shared_kmer(genotype_mapping,kmer_profile,thresh=0.95):
    genotype_kmer_pool = {}
    genotype_counts = {}
    for sample_id in kmer_profile:
        genotype = genotype_mapping[sample_id]
        if not genotype in genotype_kmer_pool:
            genotype_kmer_pool[genotype] = {}
            genotype_counts[genotype] = 0
        genotype_counts[genotype]+=1
        for kmer_id in kmer_profile[sample_id]:
            if kmer_id not in genotype_kmer_pool[genotype]:
                genotype_kmer_pool[genotype] = 0

    shared_kmers = {}
    for genotype in genotype_kmer_pool:
        count = genotype_counts[genotype]
        for kmer_id in genotype_kmer_pool[genotype]:
            value = genotype_kmer_pool[genotype][kmer_id] / count
            if value > thresh:
                if not genotype in shared_kmers:
                    shared_kmers[genotype] = []
                shared_kmers[genotype].append(kmer_id)

    return shared_kmers

def identify_diagnostic_kmer(genotype_mapping,kmer_profile,thresh=0.95):
    shared_kmers = identify_shared_kmer(genotype_mapping,kmer_profile,thresh)
    diagnostic_kmers = {}
    genotype_kmer_pool = {}
    report = {}
    for sample_id in kmer_profile:
        genotype = genotype_mapping[sample_id]
        if not genotype in genotype_kmer_pool:
            genotype_kmer_pool[genotype] = []
        genotype_kmer_pool[genotype].extend(kmer_profile[sample_id])


    for genotype_1 in shared_kmers:
        kmers = set(shared_kmers[genotype_1])
        for genotype_2 in genotype_kmer_pool:
            if genotype_1 != genotype_2:
                kmers = kmers - genotype_kmer_pool[genotype_2]
        diagnostic_kmers[genotype_1] = kmers

    return diagnostic_kmers

def qa_genotypes(genotype_mapping,kmer_profile,shared_kmers,diagnostic_kmers):
    report = {}
    unique_profiles_samples = {}
    unique_profiles_genotypes = {}
    genotypes = []
    for sample_id in kmer_profile:
        md5 = calc_md5(kmer_profile[sample_id])
        genotype = genotype_mapping[sample_id]
        genotypes.append(genotype)
        if not md5 in unique_profiles_samples:
            unique_profiles_samples[md5] = []
            unique_profiles_genotypes[md5] = []
        unique_profiles_samples[md5].append(sample_id)
        unique_profiles_genotypes[md5].append(genotype)
    genotypes = list(set(genotypes))
    for genotype in genotypes :
        report[genotype] = {
            'name': genotype,
            'num_shared_kmers': 0,
            'num_diagnostic_kmers': 0,
            'num_conflicting_profiles': 0,
            'shared_kmers': '',
            'diagnostic_kmers': '',
            'conflicting_profiles': '',
            'qc_message':''
        }

        if genotype in shared_kmers:
            report[genotype]['num_shared_kmers'] = len(shared_kmers[genotype])
            report[genotype]['shared_kmers'] = ', '.join(shared_kmers[genotype])

        if genotype in shared_kmers:
            report[genotype]['num_shared_kmers'] = len(diagnostic_kmers[genotype])
            report[genotype]['shared_kmers'] = ', '.join(diagnostic_kmers[genotype])

    genotype_conflict_counts = {}
    problem_profiles = {}
    for md5 in unique_profiles_genotypes:
        genotype_counts = {k: v for k, v in
                                  sorted(Counter(unique_profiles_genotypes[md5]).items(), key=lambda item: item[1], reverse=True)}
        num_genotypes = len(genotype_counts)
        if num_genotypes == 1:
            continue
        problem_profiles[md5] = genotype_counts
        for genotype in genotype_counts:
            #skip the conflict for the genotype if it has an exclusive kmer
            if genotype in diagnostic_kmers:
                continue

            if genotype not in genotype_conflict_counts :
                genotype_conflict_counts[genotype] = {'count':0,'profiles':[]}
            genotype_conflict_counts[genotype]['count']+=1
            genotype_conflict_counts[genotype]['profiles'].append(md5)

    for genotype in genotype_conflict_counts:
        report[genotype]['num_conflicting_profiles'] =  len(genotype_conflict_counts[genotype]['profiles'])
        report[genotype]['conflicting_profiles'] = genotype_conflict_counts[genotype]['profiles']

    for genotype in report:
        qc_message = []
        if report[genotype]['num_shared_kmers'] == 0:
            qc_message.append( 'FAIL: No shared kmers' )
        if report[genotype]['diagnostic_kmers'] == 0:
            qc_message.append( 'WARNING: No diagnostic kmers' )
        if report[genotype]['num_conflicting_profiles'] > 0:
            qc_message.append( 'WARNING: Genotype involved in conflicts' )

    return '; '.join(qc_message )

def run():
    #get arguments
    cmd_args = parse_args()
    logger = init_console_logger(2)

    #input parameters
    input_msa = cmd_args.input_msa
    input_meta = cmd_args.input_meta
    prefix = cmd_args.prefix
    ref_id = cmd_args.ref_id
    ref_gbk = cmd_args.ref_gbk
    outdir = cmd_args.outdir
    min_len = cmd_args.min_len
    max_len = cmd_args.max_len
    max_ambig = cmd_args.max_ambig
    min_complexity = cmd_args.min_complexity

    # initialize analysis directory
    if not os.path.isdir(outdir):
        logger.info("Creating analysis results directory {}".format(outdir))
        os.mkdir(outdir, 0o755)
    else:
        logger.info("Results directory {} already exits, will overwrite any results files here".format(outdir))


    #output files
    scheme_file = os.path.join(outdir,"{}-scheme.txt".format(prefix))
    genotypes_file = os.path.join(outdir,"{}-mutations.txt".format(prefix))
    genotype_dendrogram = os.path.join(outdir, "{}-dendropgram.png".format(prefix))
    reference_fasta_file = os.path.join(outdir, "{}-{}.fasta".format(prefix,ref_id))
    biohansel_fasta_file = os.path.join(outdir,"{}-biohansel-scheme.fasta".format(prefix))

    #Get the Gene features from the reference sequence
    ref_features = parse_reference_sequence(ref_gbk)

    #Read the input MSA and calculate the consensus sequence
    input_alignment = read_fasta(input_msa)
    consensus_bases = calc_consensus(input_alignment)
    consensus_seq = generate_consensus_seq(consensus_bases)

    #Identify variable positions within the alignment
    snp_positions = find_snp_positions(consensus_seq)

    sequence_deletions = {}
    for seq_id in input_alignment:
        sequence_deletions[seq_id] = find_internal_gaps(input_alignment[seq_id])

    #Identify kmers which are compatible with the user specifications around each mutation
    scheme = find_snp_kmers(input_alignment,snp_positions,consensus_bases,consensus_seq,ref_features,ref_id, \
                            min_len=min_len,max_len=max_len,max_ambig=max_ambig)
    scheme.update(find_indel_kmers(input_alignment, sequence_deletions, consensus_seq, ref_features, ref_id, \
                            min_len=min_len,max_len=max_len,max_ambig=max_ambig))

    #Identify problems with selected kmers so that the user can filter them
    scheme = qa_scheme(scheme,missing_targets=[], multiplicity_targets=[], min_len=min_len,max_len=max_len,max_ambig=max_ambig,min_complexity=min_complexity)

    #write the initial scheme to a file
    print_scheme(scheme, scheme_file)

    # write scheme to file to run biohansel on
    logger.info("Writing biohansel compatible kmer scheme to {}".format(biohansel_fasta_file))
    scheme_df = read_tsv(scheme_file)
    scheme_df = scheme_df[scheme_df['is_valid'] == True]
    scheme_to_biohansel_fasta(scheme_df,biohansel_fasta_file )
    scheme_kmer_target_info = init_kmer_targets(scheme_df)
    scheme_kmer_target_keys = list(scheme_kmer_target_info.keys())
    scheme_kmer_target_keys.sort()

    #write reference sequence to fasta
    seq = input_alignment[ref_id].replace('-','')
    fh = open(reference_fasta_file,'w')
    fh.write(">{}\n{}\n".format(ref_id,seq))


    # run kmer detection on reference sequence
    kmer_file = os.path.join(outdir, "bh.kmer.txt")
    summary_file = os.path.join(outdir, "bh.summary.txt")
    simple_file = os.path.join(outdir, "bh.simple.txt")
    logger.info("Identifying occurance of scheme kmers which are found in the reference {}".format(ref_id))
    (stdout,stderr) = bio_hansel.run_biohansel_single(biohansel_fasta_file, reference_fasta_file, kmer_file, summary_file, simple_file)
    ref_kmer_biohansel_df = read_tsv(kmer_file)
    missing_targets = list(set(scheme_kmer_target_keys) - set(ref_kmer_biohansel_df['refposition'].tolist()))
    target_counts = ref_kmer_biohansel_df['refposition'].value_counts().to_dict()
    multiplicity_targets = []
    for target in target_counts:
        if target_counts[target] > 1:
            multiplicity_targets.append(target)

    #Identify problems with selected kmers so that the user can filter them
    scheme = qa_scheme(scheme,missing_targets, multiplicity_targets, min_len=min_len,max_len=max_len,max_ambig=max_ambig,min_complexity=min_complexity)

    #write the final scheme to a file
    print_scheme(scheme, scheme_file)

    #read the metadata associations
    metadata_df = read_tsv(input_meta)
    genotype_mapping = get_genotype_mapping(metadata_df)

    #get the kmer profile for each sample
    kmer_profile = build_kmer_profiles(scheme)

    #create a plot of sample similarity for a multi-sample run
    if len(kmer_profile ) > 1:
        dist_matrix = profile_pairwise_distmatrix(kmer_profile)
        d = dendrogram_visualization()
        d.build_tree_from_dist_matrix(list(kmer_profile.keys()),dist_matrix ,genotype_dendrogram)

    #identify genotype shared kmers
    shared_kmers = identify_shared_kmer(genotype_mapping,kmer_profile,thresh=0.95)

    #identify genotype diagnostic kmers
    diagnostic_kmers = identify_diagnostic_kmer(genotype_mapping, kmer_profile, thresh=0.95)

    #Analyse the genotypes for conflicts
    genotype_report = qa_genotypes(genotype_mapping,kmer_profile,shared_kmers,diagnostic_kmers)

    # write the genotype file
    print_scheme(genotype_report, genotype_report)

