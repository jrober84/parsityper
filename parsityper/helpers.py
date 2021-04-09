from Bio import SeqIO
import pandas as pd
import logging, os, copy
from scipy.spatial.distance import cdist
from parsityper.constants import HTML_TEMPLATE_FILE, LOG_FORMAT, TYPING_SCHEMES

def init_console_logger(lvl=2):
    '''
    Controlls the level of messaging provided to the user
    :param lvl: integer indicating level
    :return: logging object
    '''
    root = logging.getLogger()

    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]

    report_lvl = logging_levels[lvl]
    root.setLevel(report_lvl)  # set root logger level

    logging.basicConfig(format=LOG_FORMAT, level=report_lvl)

    return logging.getLogger(__name__)



def read_scheme_fasta(fasta_file):
    '''
    :param fasta_file: Fasta file
    :return: dictionary of sequences indexed by id
    '''
    reference = {}
    for seq_record in SeqIO.parse(fasta_file,format='fasta'):
        reference[str(seq_record.id)] = list(str(seq_record.seq))
    return reference

def report_template():
    '''
    Reads template file into string
    :return: HTML template string
    '''
    return open(HTML_TEMPLATE_FILE, 'r').read()

def read_tsv(tsv_file):
    '''
    Reads TSV file into a pandas data frame
    :param tsv_file: tsv file path
    :return: pd dataframe
    '''
    return pd.read_csv(tsv_file,sep="\t",header=0,low_memory=False)

def scheme_to_biohansel_fasta(df,fasta_file):
    '''
    Writes a compatible kmer df into a biohansel compatible scheme fasta
    :param df: kmer scheme data frame
    :param fasta_file: fasta file path
    :return:
    '''
    out_string = []
    for row in df.itertuples():
        id = row.key
        pos_kmer = row.positive
        neg_kmer = row.negative
        out_string.append(">{}-{}\n{}\n>negative{}-{}\n{}".format(id,id,pos_kmer,id,id,neg_kmer))
    fh = open(fasta_file,'w')
    fh.write("\n".join(out_string))
    fh.close

def filter_biohansel_kmer_df(df,min_freq,min_frac):
    df = df[df['freq'] > min_freq]
    df = df[df['pos_ratio'] > min_freq]
    return df


def validate_args(cmd_args,logger):
    '''
    Validates command line parameters
    :param cmd_args: args object
    :param logger: logging object
    :return: True on success, false if any parameter fails its check
    '''
    is_valid = True
    min_cov = cmd_args.min_cov
    if min_cov < 1:
        logger.error("Error you need to specify an integer >= 1 for min_cov")
        is_valid = False

    min_cov_frac = cmd_args.min_cov_frac
    if min_cov_frac <= 0 or min_cov_frac > 1:
        logger.error("Error you need to specify a float value >0 and <= 1")
        is_valid = False

    scheme= cmd_args.scheme
    if not scheme in TYPING_SCHEMES and not os.path.isfile(scheme):
        logger.error("Error specified scheme name or file does not exist".format(scheme))
        is_valid = False

    mode = cmd_args.mode
    if mode != 'single' and mode != 'batch':
        logger.error("Error specified operating mode is invalid, enter 'single' or 'batch'")
        is_valid = False

    type = cmd_args.mode
    if type != 'single' and type != 'multi':
        logger.error("Error specified sample type is invalid, enter 'single' or 'multi'")
        is_valid = False

    R1 = cmd_args.R1
    R2 = cmd_args.R2
    SE = cmd_args.se

    if R1 is not None and R2 is not None:
        if not os.path.isfile(R1):
            logger.error("Error {} is not found".format(R1))
            is_valid = False
        if not os.path.isfile(R2):
            logger.error("Error {} is not found".format(R2))
            is_valid = False

    if SE is not None:
        if not os.path.isfile(SE):
            logger.error("Error {} is not found".format(SE))
            is_valid = False

    if (R1 is not None or R2 is not None) and SE is not None:
        logger.error("Error you have specified both paired-end reads and single-end, please specify only one")
        is_valid = False

    data_dir = cmd_args.data_dir
    if data_dir is not None:
        if not os.path.isdir(data_dir):
            logger.error("Error {} directory is not found".format(data_dir))
            is_valid = False

    if ((R1 is not None or R2 is not None) or SE is not None) and data_dir is not None:
        logger.error("Error you have specified readsets and a data directory, please specify only one")
        is_valid = False

    if (R1 is None or R2 is None) and SE is None and data_dir is None:
        logger.error("Error you need to specify either read sets or a data directory to process")
        is_valid = False

    return is_valid


def init_kmer_targets(scheme_df):
    '''

    :param scheme_df: K-mer typing scheme pandas dataframe
    :return: dict of kmer targets
    '''
    targets = {}
    fields = list(set(scheme_df.columns.tolist()) - set(['key']))
    for row in scheme_df.itertuples():
        key = row.key
        if not key in targets:
            targets[key] = {}
        for field in fields:
            targets[key][field] = getattr(row, field)

    return targets


def generate_biohansel_kmer_names(value_list):
    '''

    :param value_list: integer list
    :return: list of kmer names in biohansel format
    '''
    names = []
    for value in value_list:
        names.append("{}-{}".format(value,value))
        names.append("negative{}-{}".format(value, value))
    return names

def get_scheme_template(elements,value):
    template = dict()
    for e in elements:
        if isinstance(value,dict) or isinstance(value,list) :
            value = copy.deepcopy(value)
        template[e] = value
    return template

def calc_profile_distances(sample_profile,genotype_profiles):
    '''
    Computes the distance between the sample profile and all of the precomputed genotypes
    :param sample_profile: pandas df of kmer profile, 0,0.5,1 or nan as values for sample
    :param genotype_profiles: pandas df of kmer profile, 0,0.5,1 or nan as values for reference to compare sample
    :return: list of distances
    '''
    return cdist(sample_profile, genotype_profiles, 'euclidean')[0]

def identify_candidate_genotypes_by_dist(sample_profile,genotype_profiles,max_dist=0.5):
    '''
    Computes the distance between the sample profile and all of the precomputed genotypes and filteres out any whihc
    are too far from the distance threshold
    :param sample_profile: pandas df of kmer profile, 0,0.5,1 or nan as values for sample
    :param genotype_profiles: pandas df of kmer profile, 0,0.5,1 or nan as values for reference to compare sample
    :return: list of distances
    '''
    dists = calc_profile_distances(sample_profile,genotype_profiles)
    samples = genotype_profiles.index.tolist()
    candidates = {}
    for i in range(0,len(dists)):
        if dists[i] > max_dist:
            continue
        candidates[samples[i]] = dists[i]
    return candidates


