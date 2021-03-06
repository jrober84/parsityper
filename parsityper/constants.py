import os

LOG_FORMAT = '%(asctime)s %(name)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'

default_database_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
BATCH_HTML_REPORT = os.path.join(default_database_dir, 'batch.report.template.html')

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

TYPER_SAMPLE_SUMMARY_HEADER_BASE = [
    'sample_id',
    'scheme',
    'analysis_date',
    'sequencing_technology',
    'file_type',
    'reported_sample_type',
    'num_reads',
    'est_genome_size',
    'num_unique_kmers',
    'num_counted_unique_kmers',
    'raw_seq_files',
    'total_reads_pre',
    'total_reads_post',
    'total_bases_pre',
    'total_bases_post',
    'read_mean_len_pre',
    'read_mean_len_post',
    'read_gc_pre',
    'read_gc_post',
    'read_insert_size_peak',
    'read_duplication_rate',
    'estimated_genome_cov',
    'total_scheme_kmers',
    'num_detected_scheme_kmers',
    'ave_scheme_kmers_freq',
    'total_scheme_mutations',
    'num_detected_scheme_mutations',
    'detected_scheme_mixed_mutations',
    'detected_sample_type',
    'compatible_genotypes',
    'detected_sample_type',
    'primary_genotype',
    'primary_genotype_frac',
    'mutation_profile',
    'mutation_profile_md5',
    'mutation_profile_phrase',
    'qc_messages'

]


FIGURE_CAPTIONS = {
    'coverage_plot_caption':'Figure: Heat map plot representing the average k-mer frequency within a sample for a given genomic location',
    'mixed_target_plot_caption':'Figure: Bar plot of the target mutations which were found to have both states sorted by the number of samples where a target is mixed',
    'missing_features_plot_caption':'Figure: Bar plot of the target mutations where no kmers were found sorted by the number of samples where a target is missing',
    'genotype_abundance_plot_caption':'Figure: Bar plot of the abundance of each genotype which was found to be compatible',
    'sample_mds_caption':'Figure: MDS plot of sample similarity using k-mer based jaccard distances',
    'sample_dendro_caption':'Figure: Dendrogram of sample similarity using k-mer based jaccard distances and corresponding distance matrix shown as a heatmap',
    'feature_dendro_caption':'Figure: Dendrogram of sample similarity using k-mer based jaccard distances. Heatmap represents alt mutation kmers and their abundances'
}

iupac_replacement = {'R': 'N', 'Y': 'N', 'M': 'N', 'K': 'N', 'S': 'N', 'W': 'N', 'H': 'N', 'B': 'N', 'V': 'N', 'D': 'N'}

SCHEME_HEADER = [
    'key',
    'mutation_key',
    'mutation_type',
    'dna_name',
    'gene',
    'gene_start',
    'gene_end',
    'cds_start',
    'cds_end',
    'aa_name',
    'aa_start',
    'aa_end',
    'is_silent',
    'is_cds',
    'is_frame_shift',
    'variant_start',
    'variant_end',
    'kmer_start',
    'kmer_end',
    'target_variant',
    'target_variant_len',
    'state',
    'ref_state',
    'alt_state',
    'kseq',
    'klen',
    'homopolymer_len',
    'kmer_entropy',
    'positive_genotypes',
    'partial_genotypes',
    'is_ambig_ok',
    'is_kmer_found',
    'is_kmer_length_ok',
    'is_kmer_unique',
    'is_valid'
]

KMER_FIELDS = {
    'key': 0,
    'mutation_key': '',
    'dna_name': '',
    'variant_start': -1,
    'variant_end': -1,
    'kmer_start': -1,
    'aa_name': '',
    'gene_name': '',
    'cds_start': -1,
    'cds_end': -1,
    'aa_start': -1,
    'aa_end': -1,
    'is_silent': False,
    'is_cds': False,
    'is_frame_shift': False,
    'target_variant': '',
    'target_variant_len': '',
    'mutation_type': '',
    'state': '',
    'kseq': '',
    'klen': 0,
    'homopolymer_len': 0,
    'ref_state': '',
    'alt_state': '',
    'entropy': -1,
    'positive_genotypes': [],
    'partial_genotypes': [],
}

