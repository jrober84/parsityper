#!/usr/bin/python
import logging
from argparse import (ArgumentParser)
import os, sys, glob
import pandas as pd
from parsityper.version import __version__
from parsityper.scheme import parseScheme, constructSchemeLookups, SCHEME_HEADER
from parsityper.helpers import init_console_logger

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Analyse MSA')
    parser.add_argument('--input1', type=str, required=False, help='TSV file 1')
    parser.add_argument('--input2', type=str, required=False,help='TSV file 2')
    parser.add_argument('--input_dir', type=str, required=False, help='input directory of files to merge')
    parser.add_argument('--type', type=str, required=True, help='profile or scheme')
    parser.add_argument('--out', type=str, required=True,help='output file')
    parser.add_argument('-V', '--version', action='version', version='%(prog)s {}'.format(__version__))
    return parser.parse_args()

def merge_profile(files,outfile):
    mergedDF = pd.DataFrame()
    for profile_file in files:
        logging.info("processing file: {}".format(profile_file))
        if len(mergedDF) == 0:
            mergedDF = pd.read_csv(profile_file, sep="\t", header=0, index_col=0)
            continue
        df = pd.read_csv(profile_file, sep="\t", header=0, index_col=0)
        mergedDF = pd.merge(mergedDF, df, left_index=True, right_index=True)
    logging.info("Writing results to : {}".format(outfile))
    mergedDF.to_csv(outfile,header=True,sep="\t")

def merge_scheme(files,outfile):
    rules = {}
    partial_rules = {}
    files.sort()
    lengths = []
    for scheme_file in files:
        logging.info("processing file: {}".format(scheme_file))
        scheme_info = constructSchemeLookups(parseScheme(scheme_file))
        lengths.append(len(scheme_info['uid_to_kseq']))
        for genotype in scheme_info['genotype_rule_sets']:
            uids = scheme_info['genotype_rule_sets'][genotype]['positive_uids']
            for uid in uids:
                kmer = scheme_info['uid_to_kseq'][uid]
                if not kmer in rules:
                    rules[kmer] = []
                rules[kmer].append(genotype)

            uids = scheme_info['genotype_rule_sets'][genotype]['partial_uids']
            for uid in uids:
                kmer = scheme_info['uid_to_kseq'][uid]
                if not kmer in partial_rules:
                    partial_rules[kmer] = []
                partial_rules[kmer].append(genotype)
    fh = open(outfile,'w')
    fh.write("{}\n".format("\t".join(SCHEME_HEADER)))
    df = pd.read_csv(scheme_file, sep="\t", header=0)
    columns = df.columns.tolist()
    num_fields = len(SCHEME_HEADER)
    logging.info("Writing results to : {}".format(outfile))
    new_uid_key = 0
    for index,row in df.iterrows():
        seq = row['unalign_kseq']
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
        if seq in rules:
            positive_genotypes = rules[seq]
        positive_genotypes = list(set([str(x) for x in positive_genotypes]))
        positive_genotypes.sort()
        entry['positive_genotypes'] = ','.join(positive_genotypes)
        partial_genotypes = []
        if seq in partial_rules:
            partial_genotypes = partial_rules[seq]
        partial_genotypes = list(set([str(x) for x in partial_genotypes]))
        partial_genotypes.sort()

        entry['seq_ids'] = ''
        entry['partial_genotypes'] = ','.join(partial_genotypes)

        record = []
        for i in range(0,num_fields):
            record.append(entry[SCHEME_HEADER[i]])
        fh.write("{}\n".format("\t".join([str(x) for x in record])))
        new_uid_key+=1
    fh.close()
    return

def run():
    #input parameters
    logger = init_console_logger(2)
    cmd_args = parse_args()
    input1 = cmd_args.input1
    input2 = cmd_args.input2
    input_dir = cmd_args.input_dir
    type = str(cmd_args.type).lower()
    outfile = cmd_args.out

    #check input paramters valid
    if input1 == None and input2 == None and input_dir ==  None:
        logging.ERROR("No input files specified")
        sys.exit()
    if (input1 != None and input2 != None) and input_dir !=  None:
        logging.ERROR("You must specify two files OR an input directory")
        sys.exit()

    if (input1 == None or input2 == None) and input_dir ==  None:
        logging.ERROR("You must specify two files OR an input directory")
        sys.exit()

    if not type in ['profile','scheme']:
        logging.ERROR("You must specify the type of file to merge as a 'profile' or 'scheme'")
        sys.exit()

    files = []
    if input1 != None:
        if not os.path.isfile(input1):
            logging.ERROR("File '{}' does not exist".format(input1))
            sys.exit()
        files.append(input1)
    if input2 != None:
        if not os.path.isfile(input2):
            logging.ERROR("File '{}' does not exist".format(input2))
            sys.exit()
        files.append(input2)

    if input_dir is not None:
        if input_dir[-1] != '/':
            input_dir = "{}/".format(input_dir)
        files.extend(glob.glob("{}**".format(input_dir), recursive=True))

    valid_files = []
    for file in files:
        if os.path.isfile(file):
            valid_files.append(file)

    if len(valid_files) < 2:
        logging.ERROR("Need at least two valid files to merge found files:{}".forma(files))
        sys.exit()

    if type == 'scheme':
        merge_scheme(valid_files, outfile)
    else:
        merge_profile(valid_files,outfile)

    return

# call main function
if __name__ == '__main__':
    run()