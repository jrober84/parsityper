import pandas as pd
from argparse import (ArgumentParser, FileType)
import sys, copy

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Generate psedo-sequences from column data')
    parser.add_argument('--input', type=str, required=True,
                        help='Column formated base counts')
    parser.add_argument('--min_count', type=int, required=False,
                        help='Minimum occuances of state to be included',default=100)
    return parser.parse_args()

def generate_consensus(df):
    maxValueIndexObj = df.idxmax(axis=1)
    return maxValueIndexObj.tolist()

def create_pseuoseqs(consensus,df,min_count):
    seqs = {
        'consensus':consensus,
        'A':copy.deepcopy(consensus),
        'T': copy.deepcopy(consensus),
        'C': copy.deepcopy(consensus),
        'G': copy.deepcopy(consensus)
    }

    for row in df.itertuples():
        pos=row.pos
        if consensus[pos] == '-' or consensus[pos] == 'N':
            continue

        if row.A >= min_count:
            seqs['A'][pos] = 'A'

        if row.T >= min_count:
            seqs['T'][pos] = 'T'
        if row.C >= min_count:
            seqs['C'][pos] = 'C'
        if row.G >= min_count:
            seqs['G'][pos] = 'G'

    return seqs




def main():
    cmd_args = parse_args()
    input_file = cmd_args.input
    min_count = cmd_args.min_count
    df = pd.read_csv(input_file,sep="\t",header=0)
    seqs = create_pseuoseqs(generate_consensus(df),df,min_count)
    for seq_id in seqs:
        print(">channel_{}\n{}\n".format(seq_id,''.join(seqs[seq_id])))



# call main function
if __name__ == '__main__':
    main()