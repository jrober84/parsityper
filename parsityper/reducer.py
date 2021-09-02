from Bio import SeqIO
from argparse import (ArgumentParser, FileType)
import os, sys, re, collections, operator, math, time,base64
import pandas as pd
import hashlib, copy
from subprocess import Popen, PIPE
from Bio import GenBank
import tempfile
from parsityper.helpers import read_fasta

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Analyse MSA')
    parser.add_argument('--input_msa', type=str, required=True, help='Aligned fasta file')
    parser.add_argument('--input_meta', type=str, required=True,
                        help='tab delimited sample genotype information')
    parser.add_argument('--ref_id', type=str, required=True,
                        help='sample_id for reference sequence to use in MSA')
    parser.add_argument('--ref_gbk', type=str, required=True,
                        help='GenBank file for reference sequences')
    parser.add_argument('--outdir', type=str, required=True,
                        help='output directory')
    parser.add_argument('--max_samples', type=int, required=False, help='Maximum number of samples per genotype',default=100)
    parser.add_argument('--min_length', type=int, required=False, help='Minimum sequence length',default=-1)
    parser.add_argument('--min_freq', type=int, required=False, help='Minimum occurances of a mutation',default=1)
    parser.add_argument('--max_ambig_perc', type=float, required=False,
                        help='Maximum percentage of ambiguous characters',default=1.0)
    parser.add_argument('--min_genotype_count', type=int, required=False,
                        help='Minimum occuances of mutation to be included', default=1)
    parser.add_argument('--seed', type=int, required=False,
                        help='Seed for random generator', default=42)
    parser.add_argument('--folder_size', type=int, required=False, help='Maximum number of samples per folder',
                        default=5000)
    parser.add_argument('--n_threads', type=int, required=False,
                        help='Num threads to use', default=1)
    parser.add_argument('--train_proportion', type=float, required=False,
                        help='Num threads to use', default=0.8)
    parser.add_argument('--select_best_seqs', required=False,
                        help='Flag to toggle selecting highest quality seqs for scheme development',action='store_true')
    return parser.parse_args()


def read_metadata(file):
    return pd.read_csv(file, header=0, sep="\t")


def read_sample_mapping(df, sample_col, genotype_col):
    metadata = {}
    for index, row in df.iterrows():
        sample_id = str(row[sample_col])
        genotype = str(row[genotype_col])
        metadata[sample_id] = genotype
    return metadata


def summarize_genotypes(sample_info):
    genotype_counts = sample_info['genotype'].value_counts().reset_index(name="count")
    genotype_counts.columns = ['genotype', 'count']
    return genotype_counts


def create_fasta_folder_structure(outdir, genotypes):
    if not os.path.isdir(outdir):
        os.mkdir(outdir, 0o755)

    for genotype in genotypes:
        path = os.path.join(outdir, str(genotype))
        if not os.path.isdir(path):
            os.mkdir(path, 0o755)


def calc_md5(seq):
    seq = str(seq).encode()
    md5 = hashlib.md5()
    md5.update(seq)
    return md5.hexdigest()


def folder_tracker(genotype_counts_df):
    tracker = {}
    for row in genotype_counts_df.itertuples():
        genotype = row.genotype
        count = row.count
        tracker[genotype] = {'total': count, 'index': 0, 'tracker': 0}
    return tracker


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

def parse_mafft(out):
    lines = out.split('\n')
    seqs = {}
    id = ''
    seq = []
    for line in lines:
        line = line.strip()
        if line == '':
            continue
        if line[0] == '>':
            if id == '':
                id = line.replace('>', '')
            else:
                seqs[id] = ''.join(seq)
                seq = []
                id = line.replace('>','')
        else:
            seq.append(line)
    seqs[id] = ''.join(seq)
    return seqs


def run_mafft(in_seq,threads):
    p = Popen(['mafft',
               '--auto', '--quiet', '--thread',"{}".format(threads),
               in_seq],
              stdin=PIPE,
              stdout=PIPE,
              stderr=PIPE)
    stdout, stderr = p.communicate()
    if isinstance(stdout, bytes):
        stdout = stdout.decode()
    return (stdout, stderr)

def map_positions(seq_1,seq_2):
    look_up = {}
    seq_2_len = len(seq_2)
    seq_2_pos = 0
    for i in range(0,len(seq_1)):
        base_1 = seq_1[i]
        if base_1 == '-':
            look_up[i] = {'ref_aln_pos':i,'mapped_aln_pos':-1}
            continue
        for k in range(seq_2_pos,seq_2_len):
            base_2 = seq_2[k]
            seq_2_pos += 1
            if base_1 == base_2:
                look_up[i] = {'ref_aln_pos': i, 'mapped_aln_pos': k}
                break

    return look_up

def get_valid_positions(global_consensus,threshold=1):
    valid_positions = []
    for i in range(0, len(global_consensus)):
        data = global_consensus[i]
        count_bases = data['A'] + data['T'] + data['C'] + data['G']
        if count_bases >= threshold:
            valid_positions.append(i)
    return valid_positions

def create_alignment(ref_lookup,ref_id,ref_seq,seqs,valid_positions,n_threads=1):
    aln = {ref_id:ref_seq.upper()}
    ref_len = len(ref_seq)
    for seq_id in seqs:
        seq = seqs[seq_id]
        aln_seq = ['-'] * ref_len
        for i in range(0,len(seq)):
            base = seq[i]
            if i not in valid_positions:
                continue
            mapped_aln_pos = ref_lookup[i]['mapped_aln_pos']

            if mapped_aln_pos == -1:
                continue
            aln_seq[mapped_aln_pos] = base
        aln[seq_id] = ''.join(aln_seq)

    return aln

def mafft_add_seq(input_ref_seq,input_msa,output,n_threads):
    fh = open(output,'w')
    p = Popen(['mafft', '--add',input_ref_seq,
               '--auto', '--quiet', '--thread',"{}".format(n_threads),
               input_msa],
              stdin=PIPE,
              stdout=fh,
              stderr=PIPE)
    stdout, stderr = p.communicate()
    fh.close()
    return (stderr)


def create_training_sets(seq_info_df,proportion,max_size,seed):
    training_sets = {}
    genotypes = seq_info_df['genotype'].unique().tolist()
    seq_info_df = seq_info_df.drop_duplicates(subset="md5")
    for genotype in genotypes:
        subset = seq_info_df[seq_info_df['genotype'] == genotype].sample(frac=proportion, replace=False, random_state=seed)
        sample_ids = subset['sample_id'].tolist()
        if len(sample_ids) > max_size:
            sample_ids = sample_ids[:max_size-1]
        training_sets[genotype] = sample_ids

    return training_sets


def run():
    cmd_args = parse_args()
    fasta_file = cmd_args.input_msa
    meta_file = cmd_args.input_meta
    outdir = cmd_args.outdir
    min_genotype_count = cmd_args.min_genotype_count
    max_samples = cmd_args.max_samples
    min_length = cmd_args.min_length
    max_ambig_perc = cmd_args.max_ambig_perc
    seed = cmd_args.seed
    train_proportion = cmd_args.train_proportion
    ref_gbk = cmd_args.ref_gbk
    folder_size = cmd_args.folder_size
    min_freq = cmd_args.min_freq
    n_threads = cmd_args.n_threads

    if not os.path.isdir(outdir):
        os.mkdir(outdir, 0o755)

    fasta_outdir = os.path.join(outdir, 'fastas')
    scheme_datadir = os.path.join(outdir, 'scheme_data')
    if not os.path.isdir(scheme_datadir):
        os.mkdir(scheme_datadir, 0o755)
    benchmark_datadir = os.path.join(outdir, 'benchmark_data')
    if not os.path.isdir(benchmark_datadir):
        os.mkdir(benchmark_datadir, 0o755)

    meta_df = read_metadata(meta_file)
    genotype_count_df = summarize_genotypes(meta_df)
    genotype_count_df['genotype'] = genotype_count_df['genotype'].astype(str)
    genotype_count_df.to_csv(os.path.join(outdir, "genotype_counts.txt"), index=False, sep="\t")
    genotype_count_df = genotype_count_df[genotype_count_df['count'] >= min_genotype_count]
    valid_genotypes = genotype_count_df['genotype'].tolist()
    sample_mapping = read_sample_mapping(meta_df, 'sample_id', 'genotype')
    create_fasta_folder_structure(fasta_outdir, valid_genotypes)
    tracker = folder_tracker(genotype_count_df)

    for genotype in tracker:
        if not os.path.isdir(os.path.join(fasta_outdir, "{}/{}".format(genotype, tracker[genotype]['index']))):
            os.mkdir(os.path.join(fasta_outdir, "{}/{}".format(genotype, tracker[genotype]['index'])), 0o755)

    align_len = 0
    genotype_consensus = {}
    erroneous_seqs_fh = open(os.path.join(outdir, 'error_seqs.txt'), 'w')
    seq_report_fh = open(os.path.join(outdir, 'seq_info.txt'), 'w')
    seq_report_fh.write("sample_id\tgenotype\tmd5\tnum_seq_bases\tambig_count\tgap_count\tstatus\n")
    global_consensus = []

    for seq_record in SeqIO.parse(fasta_file, format='fasta'):
        seq = str(seq_record.seq).upper()
        length = len(seq)
        seq = re.sub(r'[^A|T|C|G|-]', 'N', seq)
        gap_count = seq.count('-')
        ambig_count = seq.count('N')
        num_seq_bases = length - (gap_count + ambig_count)
        md5 = calc_md5(seq)

        if align_len == 0:
            align_len = length
        if len(global_consensus) == 0:
            for i in range(0, length):
                global_consensus.append({'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0, '-': 0})

        id = str(seq_record.id)

        if align_len != length:
            erroneous_seqs_fh.write("{}\t{}\t{}\n".format(seq_record, length, "mismatched length"))
            continue
        elif len(id) < 1:
            erroneous_seqs_fh.write("{}\t{}\t{}\n".format(seq_record, length, "ID field malformed"))
            continue
        elif not id in sample_mapping:
            erroneous_seqs_fh.write("{}\t{}\t{}\n".format(seq_record, length, "Sample ID not in genotype associations"))
            continue

        sample_id = id
        genotype = sample_mapping[sample_id]
        ambig_perc = ambig_count / num_seq_bases
        status = 'Pass'
        if ambig_perc > max_ambig_perc or num_seq_bases < min_length or genotype not in valid_genotypes:
            status = 'Fail'

        if genotype not in genotype_consensus and status == 'Pass':
            genotype_consensus[genotype] = []
            for i in range(0, length):
                genotype_consensus[genotype].append({'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0, '-': 0})

        if status == 'Pass':
            for i in range(0, length):
                base = seq[i]
                genotype_consensus[genotype][i][base] += 1
                global_consensus[i][base] += 1
            tracker[genotype]['tracker']+=1
            if tracker[genotype]['tracker'] >= folder_size:
                tracker[genotype]['tracker'] = 0
                tracker[genotype]['index']+=1
                if not os.path.isdir(os.path.join(fasta_outdir, "{}/{}".format(genotype, tracker[genotype]['index']))):
                    os.mkdir(os.path.join(fasta_outdir, "{}/{}".format(genotype, tracker[genotype]['index'])), 0o755)

            fasta_out = os.path.join(fasta_outdir, genotype)
            fasta_out = os.path.join(fasta_out, str(tracker[genotype]['index']))
            seq_fh = open(os.path.join(fasta_out, "{}.fasta".format(sample_id)), 'w')
            seq_fh.write(">{}\n{}".format(sample_id, re.sub(r'-', '', seq)))
            seq_fh.close()


        seq_report_fh.write(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sample_id, genotype, md5, num_seq_bases, ambig_count, gap_count,status))
    erroneous_seqs_fh.close()
    seq_report_fh.close()
    seq_info_df = read_metadata(os.path.join(outdir, 'seq_info.txt'))
    seq_info_df = seq_info_df[seq_info_df['status'] == 'Pass']

    consensus_outdir = os.path.join(outdir, 'consensus')

    if not os.path.isdir(consensus_outdir):
        os.mkdir(consensus_outdir, 0o755)

    consensus_df = pd.DataFrame(global_consensus)

    consensus_df.to_csv(
        os.path.join(consensus_outdir, "{}_consensus_report.txt".format('global')))
    global_consensus_seq = generate_consensus(consensus_df)
    consensus_df['pos'] = consensus_df.index

    pseudo_seqs = create_pseuoseqs(global_consensus_seq, consensus_df, min_freq)
    cons_seqs = {}


    for genotype in genotype_consensus:
        pd.DataFrame(genotype_consensus[genotype]).to_csv(
            os.path.join(consensus_outdir, "{}_consensus_report.txt".format(genotype)))
        cons_seq = []
        for pos in genotype_consensus[genotype]:
            cons_seq.append(max(pos.items(), key=operator.itemgetter(1))[0])
        cons_seqs[genotype] = cons_seq

    with open(ref_gbk) as handle:
        for record in GenBank.parse(handle):
            gb_accession = record.accession[0]
            gb_accession_version = record.version.split('.')[1]
            genome_seq = repr(record.sequence).replace("\'",'')

    unaligned = ">{}\n{}\n>{}\n{}\n".format(gb_accession,genome_seq,"consensus",''.join(pseudo_seqs["consensus"]))
    unaligned_file = os.path.join(outdir,"unaligned.fas")

    fh = open(unaligned_file,'w')
    fh.write(unaligned)
    fh.close()
    (stdout,stderr) = run_mafft(unaligned_file,n_threads)
    os.remove(unaligned_file)
    aligned_seq = parse_mafft(stdout)
    ref_lookup = map_positions(global_consensus_seq, list(aligned_seq['consensus'].upper()))
    valid_positions = get_valid_positions(global_consensus,min_freq)


    training_sets = create_training_sets(seq_info_df, train_proportion, max_samples, seed)
    training_samples = []
    for genotype in training_sets:
        training_samples.extend(training_sets[genotype])
    test_sets = create_training_sets(seq_info_df[~seq_info_df['sample_id'].isin(training_samples)], 1-train_proportion, max_samples, seed)

    test_samples = []
    for genotype in test_sets:
        test_samples.extend(test_sets[genotype])

    train_seqs = {}
    for seq_record in SeqIO.parse(fasta_file, format='fasta'):
        seq = str(seq_record.seq).upper()
        seq = re.sub(r'[^A|T|C|G|-]', 'N', seq)
        id = str(seq_record.id)
        if id in training_samples:
            path_prefix = scheme_datadir
        elif id in test_samples:
            path_prefix = benchmark_datadir
        else:
            continue
        fasta_fh = open(os.path.join(path_prefix,"{}.fasta".format(id)), 'w')
        fasta_fh.write(">{}\n{}\n".format(id,seq.replace('-','')))
        fasta_fh.close()
        train_seqs[id] = seq

    pseudo_seqs.update(cons_seqs)
    pseudo_seqs.update(train_seqs)
    ref_seq_path = os.path.join(outdir, "ref.unaligned.fasta")
    input_msa_path = os.path.join(outdir, "pseudo.unaligned.fasta")
    fasta_fh = open(input_msa_path, 'w')
    for id in pseudo_seqs:
        if isinstance(pseudo_seqs[id],list):
            pseudo_seqs[id] = ','.join(pseudo_seqs[id])
        fasta_fh.write(">{}\n{}\n".format(id, pseudo_seqs[id]))
    fasta_fh.close()
    fasta_fh = open(ref_seq_path, 'w')
    fasta_fh.write(">{}\n{}\n".format(gb_accession, genome_seq))
    output_msa = os.path.join(consensus_outdir, "allgenotype_consensus.fasta")
    mafft_add_seq(ref_seq_path, input_msa_path, output_msa, n_threads)

# call main function
if __name__ == '__main__':
    run()
