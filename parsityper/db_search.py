import pandas as pd
from parsityper.helpers import read_tsv

def process_strain_db(tsv_file):
    data = {}
    df = read_tsv(tsv_file)
    for row in df.itertuples():
        sample_id = row.sample_id
        genotype = row.genotype
        md5 = row.kmer_st_md5
        profile = row.kmer_profile_st.replace(' ', '').split(',')
        phrase = row.kmer_phrase
        if not md5 in data:
            data[md5] = {
                'samples':{},
                'profile':profile,
                'phrase':phrase
            }
        data[md5]['samples'][sample_id] = genotype
    return data

def get_neighbours(genotype_data,kmer_profile_st):
