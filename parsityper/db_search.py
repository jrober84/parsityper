import pandas as pd
from parsityper.helpers import read_tsv

def process_strain_db(tsv_file):
    data = {}
    df = read_tsv(tsv_file)
    for row in df.itertuples():
        sample_id = str(row.sample_id)
        genotype = str(row.genotype)
        md5 = row.md5
        if isinstance(row.profile,float):
            profile = []
        else:
            profile = row.profile.replace(' ', '').split(';')
        profile = [str(x) for x in profile]
        phrase = row.phrase
        if not md5 in data:
            data[md5] = {
                'samples':{},
                'profile':set(profile),
                'phrase':phrase,
                'genotype':genotype
            }
        data[md5]['samples'][sample_id] = genotype
    return data

def get_neighbours(profile_st,profile_md5,kmer_profiles,threshold=0.05):
    profile = set(profile_st.replace(' ', '').split(';'))
    neighbors = []
    for md5 in kmer_profiles:
        target = kmer_profiles[md5]['profile']
        if md5 == profile_md5:
            jaccard = 0
        elif len(profile) > 0 or len(target) > 0:
            jaccard = 1 - (len(profile & target ) / len(profile | target))
        else:
            jaccard = 1
        if jaccard < threshold:
            neighbors.append(md5)
    return neighbors

def convert_st_to_profile(profile_st,scheme_targets):
    profile = []
    for i in range(0,len(scheme_targets)):
        if scheme_targets[i] in profile_st:
            state = 1
            if "-{}".format(scheme_targets[i]) in profile_st:
                state = 0.5
        else:
            state = 0
        profile.append(state)
    return profile