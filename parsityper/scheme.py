import logging
import sys,time

import pandas as pd
from parsityper.constants import  SCHEME_HEADER


def parseScheme(scheme_file):
    scheme = {}
    df = pd.read_csv(scheme_file,sep="\t",header=0,low_memory=False)
    uid = -1
    for row in df.itertuples():
        mutation_key = row.mutation_key
        dna_name = row.dna_name
        positive_genotypes = row.positive_genotypes
        if isinstance(positive_genotypes,float):
            positive_genotypes = []
        else:
            positive_genotypes = positive_genotypes.split(',')

        partial_genotypes = row.partial_genotypes
        if isinstance(partial_genotypes,float):
            partial_genotypes = []
        else:
            partial_genotypes = partial_genotypes.split(',')
        state = row.state
        seq = row.kseq
        gene = row.gene_name
        variant_start = row.variant_start
        variant_end = row.variant_end
        ref_state = row.ref_state
        alt_state = row.alt_state
        entropy = row.entropy

        if uid !=-1:
            if row.key - uid != 1:
                logging.error("The scheme unique id column is not sequential offending keys {}..{}".format(uid,row.key))
                sys.exit()
        uid = row.key
        aa_name = row.aa_name
        is_cds = row.is_cds
        if not mutation_key in scheme:
            scheme[mutation_key] = {
                'ref':{},
                'alt':{}
            }

        scheme[mutation_key][state][uid] = row._asdict()

    #clear nan values
    for mutation_key in scheme:
        for state in scheme[mutation_key]:
            for uid in scheme[mutation_key][state]:
                for field in scheme[mutation_key][state][uid]:
                    value = str(scheme[mutation_key][state][uid][field])
                    if value == 'nan':
                        scheme[mutation_key][state][uid][field] = ''
    return scheme

def constructSchemeLookups(scheme):

    profiles = {
        'max_variant_positions':0,
        'gene_features':[],
        'genotypes':[],
        'genotype_rule_sets':{},
        'kmer_to_genotypes':{},
        'uid_to_state':{},
        'inf_alt_uids':set(),
        'inf_ref_uids':set(),
        'uid_to_kseq':{},
        'kseq_to_uids':{},
        'uid_to_mutation':{},
        'uid_to_dna_name': {},
        'uid_to_aa_name': {},
        'uid_to_gene_feature':{},
        'uid_to_entropy': {},
        'mutation_to_uid': {},
        'kmer_profiles':{},
        'mutation_profiles': {},
        'mutation_positions':{},
        'min_kmer_len':100000,
        'max_kmer_len':0
    }
    kmers = {}
    kmer_profiles = {}
    mutation_profiles = {}
    mutations = list(scheme.keys())

    #traverse to get genotypes and kmers
    for mutation_key in scheme:
        profiles['mutation_to_uid'][mutation_key] = []
        for state in scheme[mutation_key]:
            for uid in scheme[mutation_key][state]:
                if scheme[mutation_key][state][uid]['is_cds']:
                    feature_name = scheme[mutation_key][state][uid]['gene_name']
                else:
                    feature_name = 'intergenic'
                variant_end = scheme[mutation_key][state][uid]['variant_end']
                if profiles['max_variant_positions'] < variant_end:
                    profiles['max_variant_positions'] = variant_end
                profiles['uid_to_entropy'][uid] = scheme[mutation_key][state][uid]['entropy']
                profiles['uid_to_gene_feature'][uid] = feature_name
                profiles['gene_features'].append(feature_name)
                profiles['gene_features'] = list(set(profiles['gene_features']))
                seq = scheme[mutation_key][state][uid]['kseq']
                kmers[uid] = seq
                klen = len(seq)
                if profiles['min_kmer_len'] > klen:
                    profiles['min_kmer_len'] = klen
                if profiles['max_kmer_len'] < klen:
                    profiles['max_kmer_len'] = klen

                if not seq in profiles['kseq_to_uids']:
                    profiles['kseq_to_uids'][seq] = []
                profiles['uid_to_state'][uid] = state
                if len(scheme[mutation_key][state][uid]['positive_genotypes']) > 0:
                    if state == 'alt':
                        profiles['inf_alt_uids'].add(uid)
                    else:
                        profiles['inf_ref_uids'].add(uid)
                profiles['uid_to_dna_name'][uid] = scheme[mutation_key][state][uid]['dna_name']
                profiles['uid_to_aa_name'][uid] = scheme[mutation_key][state][uid]['aa_name']
                profiles['kseq_to_uids'][seq].append(uid)
                profiles['uid_to_mutation'][uid] = mutation_key
                profiles['mutation_to_uid'][mutation_key].append(uid)
                genotypes = scheme[mutation_key][state][uid]['positive_genotypes'].split(',') + scheme[mutation_key][state][uid]['partial_genotypes'].split(',')
                profiles['kmer_to_genotypes'][uid] = genotypes
                for g in genotypes:
                    if len(g) == 0:
                        continue
                    kmer_profiles[g] = []
                    mutation_profiles[g] = []

    genotypes = list(kmer_profiles.keys())

    #init the kmer profiles
    for g in kmer_profiles:
        kmer_profiles[g] = [-1] * len(kmers)
        mutation_profiles[g] = [-1] * len(mutations)
        profiles['genotype_rule_sets'][g] = {'positive_uids':[],'positive_ref':[],'positive_alt':[],'partial_uids':[],'partial_ref':[],'partial_alt':[]}


    #populate the profiles
    for i in range(0,len(mutations)):
        mutation_key = mutations[i]
        profiles['mutation_positions'][mutation_key] = i
        for state in scheme[mutation_key]:
            for uid in scheme[mutation_key][state]:
                pos = scheme[mutation_key][state][uid]['positive_genotypes'].split(',')
                par = scheme[mutation_key][state][uid]['partial_genotypes'].split(',')
                for g in pos:
                    if g == '':
                        continue
                    profiles['genotype_rule_sets'][g]['positive_uids'].append(uid)
                    if state == 'ref':
                        mutation_profiles[g][i] = 0.0
                        profiles['genotype_rule_sets'][g]['positive_ref'].append(uid)
                    else:
                        profiles['genotype_rule_sets'][g]['positive_alt'].append(uid)
                        mutation_profiles[g][i] = 1.0
                    kmer_profiles[g][uid] = 1.0
                for g in par:
                    if g == '':
                        continue
                    profiles['genotype_rule_sets'][g]['partial_uids'].append(uid)
                    if state == 'ref':
                        profiles['genotype_rule_sets'][g]['partial_ref'].append(uid)
                    else:
                        profiles['genotype_rule_sets'][g]['partial_alt'].append(uid)

                    mutation_profiles[g][i] = 0.5
                    kmer_profiles[g][uid] = 0.5

    profiles['genotypes'] = genotypes
    profiles['uid_to_kseq'] = kmers
    profiles['kmer_profiles'] = kmer_profiles
    profiles['mutation_profiles'] = mutation_profiles

    return profiles

def detectAmbigGenotypes(scheme_info):
    genotypes = scheme_info['genotypes']
    kmer_profiles = scheme_info['kmer_profiles']
    mutation_profiles = scheme_info['mutation_profiles']
    conflict_genotype_profiles = {'kmers':{},'mutations':{}}
    for i in range(0,len(genotypes)):
        stime = time.time()
        geno_1 = genotypes[i]
        kProfile_1 =  kmer_profiles[geno_1]
        mProfile_1 = mutation_profiles[geno_1]
        for k in range(i+1,len(genotypes)):
            geno_2 = genotypes[k]
            kProfile_2 = kmer_profiles[geno_2]
            mProfile_2 = mutation_profiles[geno_2]
            kDiff = 0
            for j in range(0,len(kProfile_1)):
                if (kProfile_1[j] == 1 and kProfile_2[j] == 0) or (kProfile_1[j] == 0 and kProfile_2[j] == 1 ):
                    kDiff+=1
            mDiff = 0
            for j in range(0,len(mProfile_1)):
                if (mProfile_1[j] == 1 and mProfile_2[j] == 0 )or (mProfile_1[j] == 0 and mProfile_2[j] == 1 ):
                    mDiff+=1
            if kDiff == 0:
                if not geno_1 in conflict_genotype_profiles['kmers']:
                    conflict_genotype_profiles['kmers'][geno_1] = []
                conflict_genotype_profiles['kmers'][geno_1].append(geno_2)
            if mDiff == 0:
                if not geno_1 in conflict_genotype_profiles['mutations']:
                    conflict_genotype_profiles['mutations'][geno_1] = []
                conflict_genotype_profiles['mutations'][geno_1].append(geno_2)
    return conflict_genotype_profiles






