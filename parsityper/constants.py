import os

LOG_FORMAT = '%(asctime)s %(name)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'

default_database_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')

RUN_INFO = {
    'Program Parameters': '',
    'Min coverage':0,
    'Min coverage fraction': 0,
    'Input data': '',
    'Analysis Date':'',
    'Analysis Start time': '',
    'Analysis End time': '',
    'Primer set name':'',
    'Number of Samples':0,
    'Average number of kmers per sample': 0,
    'Stdev number of kmers per sample': 0,
    'Average kmer coverage per sample':0,
    'Stdev kmer coverage per sample':0,
    'Average mixed kmers per sample': 0,
    'Stdev mixed kmers per sample': 0,
    'Number of QC Pass samples':0,
    'Number of QC Warning samples':0,
    'Number of QC Fail samples':0,
    'Number of unique genotypes':0,
    'Number of distinct profiles':0,
}


PRIMER_SCHEMES = {
    'arctic_v1':os.path.join(default_database_dir,'arctic_v1_primers.txt'),
    'arctic_v2':os.path.join(default_database_dir,'arctic_v2_primers.txt'),
    'arctic_v3':os.path.join(default_database_dir,'arctic_v3_primers.txt'),
    'freed_v1':os.path.join(default_database_dir,'freed_v1_primers.txt'),
    'resende_v1':os.path.join(default_database_dir,'resende_v1_primers.txt')
}

TYPING_SCHEMES = {
    'SARS-COV-2_v1': os.path.join(default_database_dir,'sars-cov-2-kmers.v1.scheme.txt'),
    'SARS-COV-2_v2': os.path.join(default_database_dir,'sars-cov-2-kmers.v2.scheme.txt'),
    'SARS-COV-2_v3': '',
    'NEXTCLADE-COV-2_v1':os.path.join(default_database_dir,'sars-cov-2-kmers.nextclade.v1.scheme.txt')

}

PANGOLIN_KMER_PROFILES = {
    'SARS-COV-2_v1': os.path.join(default_database_dir,'sars-cov-2-kmers.v1.genotype_profiles.txt'),
    'SARS-COV-2_v2': '',
    'SARS-COV-2_v3': '',
}

HTML_TEMPLATE_FILE = os.path.join(default_database_dir,'report.template.html')

bases_dict = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T'], }


IUPAC_LOOK_UP = {
        'A': 'A',
        'T': 'T',
        'C': 'C',
        'G': 'G',
        'AT': 'W',
        'AC': 'M',
        'AG': 'R',
        'CT': 'Y',
        'GT': 'K',
        'CG': 'S',
        'CGT': 'B',
        'AGT': 'D',
        'ACT': 'H',
        'ACG': 'V',
        'ACGT': 'N'
    }

NEGATE_BASE_IUPAC = {'A':'B', 'C':'D','G':'H','T':'V'}


KMER_DF_HEADERS = [
    'sample',
    'kmername',
    'seq',
    'freq',
    'match_index',
    'is_pos_kmer',
    'file_path',
]


SCHEME_HEADER = [
'key',
'index',
'group_id',
'dna_name',
'aa_name',
'type',
'gene',
'gene_start',
'gene_end',
'cds_start',
'cds_end',
'is_silent',
'is_frame_shift',
'ref_aa',
'alt_aa',
'variant_start',
'variant_end',
'kmer_start',
'kmer_end',
'variant_state',
'ref_state',
'is_positive',
'positive_seqs',
'partial_positive_seqs',
'kmer_len',
'count_ambig',
'gc',
'kmer_seq',
'is_length_ok',
'is_ambig_ok',
'is_unique',
'is_detectable',
'is_valid'
]
