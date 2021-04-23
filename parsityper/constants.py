import os

LOG_FORMAT = '%(asctime)s %(name)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'

default_database_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')


PRIMER_SCHEMES = {
    'arctic_v1':os.path.join(default_database_dir,'arctic_v1_primers.txt'),
    'arctic_v2':os.path.join(default_database_dir,'arctic_v2_primers.txt'),
    'arctic_v3':os.path.join(default_database_dir,'arctic_v3_primers.txt'),
    'freed_v1':os.path.join(default_database_dir,'freed_v1_primers.txt')
}

TYPING_SCHEMES = {
    'SARS-COV-2_v1': os.path.join(default_database_dir,'sars-cov-2-kmers.v1.scheme.txt'),
    'SARS-COV-2_v2': '',
    'SARS-COV-2_v3': '',
}

PANGOLIN_KMER_PROFILES = {
    'SARS-COV-2_v1': os.path.join(default_database_dir,'sars-cov-2-kmers.v1.genotype_profiles.txt'),
    'SARS-COV-2_v2': '',
    'SARS-COV-2_v3': '',
}

HTML_TEMPLATE_FILE = os.path.join(default_database_dir,'report.template.html')