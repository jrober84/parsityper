import sys,time

import pandas as pd

SCHEME_HEADER = [
                'key',
                'mutation_key',
                'index',
                'dna_name',
                'align_variant_start',
                'align_variant_end',
                'unalign_variant_start',
                'unalign_variant_end',
                'align_kmer_start',
                'align_kmer_end',
                'unalign_kmer_start',
                'unalign_kmer_end',
                'target_variant',
                'target_variant_len',
                'mutation_type',
                'state',
                'unalign_kseq',
                'unalign_klen',
                'complexity',
                'gene',
                'gene_start',
                'gene_end',
                'cds_start',
                'cds_end',
                'aa_name',
                'aa_start',
                'aa_end',
                'ref_state',
                'alt_state',
                'is_silent',
                'is_cds',
                'is_frame_shift',
                'seq_ids',
                'positive_genotypes',
                'partial_genotypes',
                'is_kmer_found',
                'is_kmer_length_ok',
                'is_kmer_unique',
                'is_kmer_complexity_ok',
                'is_valid'
]

def parseScheme(scheme_file):
    scheme = {}
    df = pd.read_csv(scheme_file,sep="\t",header=0)
    for row in df.itertuples():
        mutation_key = row.mutation_key
        dna_name = row.dna_name
       # aa_name = row.aa_name
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
        seq = row.unalign_kseq
        gene = row.gene
        variant_start = row.unalign_variant_start
        variant_end = row.unalign_variant_end
        ref_state = row.ref_state
        alt_state = row.alt_state
        seq_ids = row.seq_ids
        uid = row.key
        if not mutation_key in scheme:
            scheme[mutation_key] = {
                'ref':{},
                'alt':{}
            }

        scheme[mutation_key][state][uid] = {
            'dna_name':dna_name,
           # 'aa_name':aa_name,
            'gene':gene,
            'variant_start':variant_start,
            'variant_end':variant_end,
            'ref_state':ref_state,
            'alt_state':alt_state,
            'seq_ids':seq_ids,
            'positive_genotypes':positive_genotypes,
            'partial_genotypes':partial_genotypes,
            'seq':seq
        }

    return scheme

def constructSchemeLookups(scheme):
    profiles = {
        'genotypes':[],
        'uid_to_state':{},
        'uid_to_kseq':{},
        'kseq_to_uids':{},
        'uid_to_mutation':{},
        'mutation_to_uid': {},
        'kmer_profiles':{},
        'mutation_profiles': {},
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
                seq = scheme[mutation_key][state][uid]['seq']
                kmers[uid] = seq
                klen = len(seq)
                if profiles['min_kmer_len'] > klen:
                    profiles['min_kmer_len'] = klen
                if profiles['max_kmer_len'] < klen:
                    profiles['max_kmer_len'] = klen

                if not seq in profiles['kseq_to_uids']:
                    profiles['kseq_to_uids'][seq] = []
                profiles['uid_to_state'][uid] = state
                profiles['kseq_to_uids'][seq].append(uid)
                profiles['uid_to_mutation'][uid] = mutation_key
                profiles['mutation_to_uid'][mutation_key].append(uid)
                genotypes = scheme[mutation_key][state][uid]['positive_genotypes'] + scheme[mutation_key][state][uid]['partial_genotypes']
                for g in genotypes:
                    kmer_profiles[g] = []
                    mutation_profiles[g] = []

    genotypes = list(kmer_profiles.keys())

    #init the kmer profiles
    for g in kmer_profiles:
        kmer_profiles[g] = [0] * len(kmers)
        mutation_profiles[g] = [0] * len(mutations)


    #populate the profiles
    for i in range(0,len(mutations)):
        mutation_key = mutations[i]
        for state in scheme[mutation_key]:
            for uid in scheme[mutation_key][state]:
                pos = scheme[mutation_key][state][uid]['positive_genotypes']
                par = scheme[mutation_key][state][uid]['partial_genotypes']

                for g in pos:
                    if state == 'ref':
                        mutation_profiles[g][i] = 0
                    else:
                        mutation_profiles[g][i] = 1
                    kmer_profiles[g][uid] = 1
                for g in par:
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
        print(time.time() - stime)
    return conflict_genotype_profiles






