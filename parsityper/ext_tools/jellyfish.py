import os.path
import pandas as pd
from parsityper.ext_tools import run_command


def run_jellyfish_count(seq_file,out_file,k=21,n_threads=1):
    cmd = "jellyfish count --text -t {} -m {} -s 100M -C -o {} {}".format(n_threads,k,out_file,seq_file)
    (stdout,stderr) = run_command(cmd)
    return (stdout,stderr)

def parse_jellyfish_counts(file):
    return pd.read_csv(file,skiprows=1, header=None,names=['kmer','count'],sep=' ')


