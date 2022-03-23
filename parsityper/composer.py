import shutil

from Bio import SeqIO
from argparse import (ArgumentParser, FileType)
import os, sys, re, collections, operator, math, time,base64
import pandas as pd
import hashlib, copy
from subprocess import Popen, PIPE
from parsityper.helpers import read_fasta
from parsityper.version import __version__
from functools import partial
from mimetypes import guess_type
import gzip
from parsityper.ext_tools.jellyfish import run_jellyfish_count, parse_jellyfish_counts
from parsityper.constants import iupac_replacement
from parsityper.kmerSearch.kmerSearch import init_automaton_dict

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Analyse MSA')
    parser.add_argument('--input_msa', type=str, required=True, help='Aligned fasta file')
    parser.add_argument('--input_meta', type=str, required=True,
                        help='tab delimited sample genotype information')
    parser.add_argument('--ref_fasta', type=str, required=True,
                        help='reference fasta')
    parser.add_argument('--mask', type=str, required=False,
                        help='TSV file of start, end regions to mask in alignment')
    parser.add_argument('--outdir', type=str, required=True,
                        help='output directory')
    parser.add_argument('--klen', type=int, required=False, help='kmer length',
                        default=18)
    parser.add_argument('--min_kmer_count', type=int, required=False, help='Mininumum frequency for kmer',
                        default=1)
    parser.add_argument('--max_homo', type=int, required=False,
                        help='Absolute maximum of homopolymer run default=20% of kmer len')
    parser.add_argument('--max_samples', type=int, required=False, help='Maximum number of samples per genotype',default=100)
    parser.add_argument('--min_length', type=int, required=False, help='Minimum sequence length',default=-1)
    parser.add_argument('--max_ambig', type=int, required=False,
                        help='Maximum number of ambiguous characters',default=0)
    parser.add_argument('-d', type=str, required=False,
                        help='fasta header delimiter', default='|')
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
    parser.add_argument('-V', '--version', action='version', version='%(prog)s {}'.format(__version__))

    return parser.parse_args()


def read_metadata(file):
    return pd.read_csv(file, header=0, sep="\t")

def init_consensus(seq):
    seq_len = len(seq)
    consensus = []
    for i in range(0, seq_len):
        consensus.append({'-': 0,'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0, })
    return consensus

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

def calc_md5(seq):
    seq = str(seq).encode()
    md5 = hashlib.md5()
    md5.update(seq)
    return md5.hexdigest()

def prep_seqs(fasta_file, genotype_mapping, delimeter, out_dir, iupac_replacement,max_ambig=0,min_length=0):
    out_file = os.path.join(out_dir,"unaligned.fasta")
    trans = str.maketrans(''.join(iupac_replacement.keys()), ''.join(iupac_replacement.values()))
    encoding = guess_type(fasta_file)[1]
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
    outFH = open(out_file,'w')
    unique_seqs = {}
    seq_report = {}
    genotypes = list(set(genotype_mapping.values()))
    genotype_consensus = {}
    genotype_files = {}
    genotype_fh = {}
    for genotype in genotypes:
        genotype_files[genotype] = os.path.join(out_dir,"{}.fasta".format(genotype))
        genotype_fh[genotype] = open(genotype_files[genotype],'w')

    with _open(fasta_file) as f:
        seq_record = next(SeqIO.parse(f, 'fasta'))
        seq = str(seq_record.seq).upper()
        consensus = init_consensus(seq)
        align_len = len(seq)
        seq_aln_range = range(0,align_len)
        for genotype in genotypes:
            genotype_consensus[genotype] = copy.deepcopy(consensus)


    with _open(fasta_file) as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            header = str(seq_record.id).split(delimeter)
            seq = str(seq_record.seq).upper()
            seq = seq.translate(trans)
            if len(seq) != align_len:
                print("Error sequences are not aligned {}".format(id))
                sys.exit()
            id = ''
            for element in header:
                if element in genotype_mapping:
                    id = element
            if id == '':
                continue
            genotype = genotype_mapping[id]
            for i in seq_aln_range:
                consensus[i][seq[i]]+=1
                genotype_consensus[genotype][i][seq[i]]+=1
            ambig_count = seq.count('N')
            seq_aln = seq
            seq = seq.replace('-','')
            is_len_ok = True
            is_ambig_ok = True
            is_duplicate = False
            if len(seq) < min_length:
                is_len_ok = False
            if ambig_count > max_ambig:
                is_ambig_ok = False
            md5 = calc_md5(seq)
            if md5 in unique_seqs:
                is_duplicate = True
            else:
                unique_seqs[md5] = id
            status = 'PASS'
            if is_duplicate or not is_ambig_ok or not is_len_ok:
                status = 'FAIL'
            seq_report[id] = {'sample_id':id,'genotype':genotype,'len':len(seq),'ambig':ambig_count,
                              'is_len_ok':is_len_ok,'is_ambig_ok':is_ambig_ok,'is_duplicate':is_duplicate,'status':status}
            if status == 'PASS':
                genotype_fh[genotype].write(">{}\n{}\n".format(id,seq_aln))
                outFH.write(">{}\n{}\n".format(id, seq))
    outFH.close()
    for fh in genotype_fh:
        genotype_fh[fh].close()
    return {'unique_samples':unique_seqs,'sample_report':seq_report,'genotype_consensus':genotype_consensus,'global_consensus':consensus}

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

def create_consensus_seq(consensus):
    seq = []
    for i in range(0,len(consensus)):
        seq.append(max(consensus[i], key=consensus[i].get))
    return ''.join(seq)

def init_mask_regions(aln_len,seq_aln_lookup,mask_file):
    aln_mask = [0] * aln_len
    df = read_metadata(mask_file)
    for index,row in df.iterrows():
        start = seq_aln_lookup[row.start]
        end = seq_aln_lookup[row.end]
        for i in range(start,end+1):
            aln_mask[i] = 1
    return aln_mask

def init_skip_regions(consensus_seqs):
    key = list(consensus_seqs.keys())[0]
    aln_len = len(consensus_seqs[key])
    aln_skip = [0] * aln_len
    aln_range = range(0,aln_len)
    for seq_id in consensus_seqs:
        for i in aln_range:
            if consensus_seqs[seq_id][i] != '-':
                aln_skip[i] = 1
    for i in aln_range:
        if aln_skip[i] == 1:
            aln_skip[i] = 0
        else:
            aln_skip[i] = 1

    return aln_skip

def create_aln_pos_from_unalign_pos(aln_seq):
    unalign_seq = aln_seq.replace('-', '')
    aln_len = len(aln_seq)
    unaln_len = len(unalign_seq)
    lookup = [-1] * unaln_len
    pos = 0
    for i in range(0, unaln_len):
        for k in range(pos, aln_len):
            if unalign_seq[i] == aln_seq[k]:
                lookup[i] = k
                pos = k + 1
                break
    return lookup

def mask_alignment_files(in_fasta_file,out_fasta_file,aln_mask,del_mask,strip=False):
    outFH = open(out_fasta_file,'w')
    encoding = guess_type(in_fasta_file)
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
    aln_len=len(aln_mask)
    with _open(in_fasta_file) as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            id = str(seq_record.id)
            seq = list(str(seq_record.seq).upper())
            seq += ['-'] * (aln_len - len(seq))
            for idx, value in enumerate(aln_mask):
                if value == 1:
                    seq[idx] = '-'

            if strip:
                seq = list(''.join(seq).replace('-',''))
            filt = []
            for idx, value in enumerate(del_mask):
                if value == 0:
                    filt.append(seq[idx])

            outFH.write(">{}\n{}\n".format(id,''.join(filt)))
    outFH.close()

def perform_aln_masking_file(in_fasta_file,out_fasta_file,ref_id,mask_file):
    encoding = guess_type(in_fasta_file)
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
    consensus = None
    ref_aln = ''
    seq_aln_range = None
    num_samples = 0
    with _open(in_fasta_file) as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            id = str(seq_record.id)
            seq = list(str(seq_record.seq).upper())
            if consensus is None:
                consensus = init_consensus(seq)
                aln_len = len(seq)
                seq_aln_range = range(0,aln_len)
            if id == ref_id:
                ref_aln = ''.join(seq)
            for i in seq_aln_range:
                consensus[i][seq[i]]+=1
            num_samples+=1
    if consensus is None:
        return False
    ref_aln_lookup = create_aln_pos_from_unalign_pos(ref_aln)
    if mask_file is not None:
        aln_mask = init_mask_regions(aln_len, ref_aln_lookup, mask_file)
    else:
        aln_mask = [0] * aln_len
    del_mask = [0] * aln_len
    for idx, value in enumerate(aln_mask):
        if ref_aln[idx] == '-':
            aln_mask[idx] = 1
            del_mask[idx] = 1
        else:
            break
    for idx, value in reversed(list(enumerate(aln_mask))):
        if ref_aln[idx] == '-':
            aln_mask[idx] = 1
            del_mask[idx] = 1
        else:
            break

    for i in seq_aln_range:
        if consensus[i]['-'] == num_samples:
            del_mask[i] = 1

    outFH = open(out_fasta_file,'w')
    #write ref as first sequence
    seq = list(ref_aln)
    filt = []
    for idx, value in enumerate(del_mask):
        if value == 0:
            filt.append(seq[idx])

    outFH.write(">{}\n{}\n".format(ref_id, ''.join(filt)))

    with _open(in_fasta_file) as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            id = str(seq_record.id)
            seq = list(str(seq_record.seq).upper())
            if id == ref_id:
                continue

            for idx, value in enumerate(aln_mask):
                if value == 1:
                    seq[idx] = '-'
            filt = []
            for idx, value in enumerate(del_mask):
                if value == 0:
                    filt.append(seq[idx])

            outFH.write(">{}\n{}\n".format(id,''.join(filt)))
    outFH.close()
    return True

def mask_alignment_dict(seqs,aln_mask,del_mask,ref_id):
    for id in seqs:
        seq = list(seqs[id])
        if id != ref_id:
            for idx, value in enumerate(aln_mask):
                if value == 1:
                    seq[idx] = '-'
        filt = []
        for idx, value in enumerate(del_mask):
            if value == 0:
                filt.append(seq[idx])
        seqs[id] = ''.join(seq)
    return seqs

def perform_aln_masking(genotypes,path_prefix,aln_mask,del_mask,strip=False):
    for genotype in genotypes:
        in_fasta_file = os.path.join(path_prefix,"{}.fasta".format(genotype))
        if strip:
            out_fasta_file = os.path.join(path_prefix,"{}.mask.unaln.fasta".format(genotype))
        else:
            out_fasta_file = os.path.join(path_prefix, "{}.mask.aln.fasta".format(genotype))
        mask_alignment_files(in_fasta_file, out_fasta_file, aln_mask, del_mask,strip)


def write_fasta_dict(seqs,out_fasta):
    outFH = open(out_fasta,'w')
    for seq_id in seqs:
        outFH.write(">{}\n{}\n".format(seq_id,seqs[seq_id]))
    outFH.close()

def get_kmer_counts(jellyfish_file):
    kmers = {}
    kmer_df = parse_jellyfish_counts(jellyfish_file)
    for row in kmer_df.itertuples():
        kmer = row.kmer
        count = row.count
        if not kmer in kmers:
            kmers[kmer] = 0
        kmers[kmer] += count
    return kmers

def calc_homopolymers(seq):
    longest = 0
    for b in ['A', 'T', 'C', 'C']:
        matches = re.findall("{}+".format(b), seq)
        for m in matches:
            length = len(m)
            if length > longest:
                longest = length
    return longest

def filter_kmers(kMers, min_count, max_count, max_ambig, max_homo):
    filtered = {}
    for kmer in kMers:
        count = kMers[kmer]
        if count < min_count or count >= max_count:
            continue
        n_count = kmer.count('N')
        if n_count > max_ambig:
            continue
        homo_len = calc_homopolymers(kmer)
        if homo_len > max_homo:
            continue
        filtered[kmer] = count
    return filtered

def get_min_seqs_per_genotype(in_fasta_file,out_fasta_file,kmer_index,aho):
    needed_kmers = set(list(kmer_index.values()))
    selected_seqs = []
    out_fh = open(out_fasta_file,'w')
    encoding = guess_type(in_fasta_file)
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
    with _open(in_fasta_file) as f:

        for seq_record in SeqIO.parse(f, 'fasta'):
            found_kmers = set()
            if len(needed_kmers) == 0:
                break
            id = str(seq_record.id)
            seq = str(seq_record.seq)
            for idx, (kIndex, kmer_seq, is_revcomp) in aho.iter(seq.replace('-','')):
                kIndex = int(kIndex)
                found_kmers.add(kmer_index[kIndex])
            if len(found_kmers) > 0:
                selected_seqs.append(id)
                needed_kmers = needed_kmers - found_kmers
                kmer_index = {}
                i = 0
                for kmer in needed_kmers:
                    kmer_index[i] = kmer
                    i += 1
                aho = init_automaton_dict(kmer_index)
                out_fh.write(">{}\n{}\n".format(id,seq))
    out_fh.close()
    return {'selected_seqs':selected_seqs,'needed_kmers':needed_kmers}

def write_individual_files(dataset,sample_mapping,test_dir,train_dir,iupac_replacement,in_fasta_file):
    trans = str.maketrans(''.join(iupac_replacement.keys()), ''.join(iupac_replacement.values()))
    encoding = guess_type(in_fasta_file)
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
    file_paths = {}
    for genotype in dataset:
        file_paths[genotype] = {
            'train':os.path.join(train_dir,"{}".format(genotype)),
            'test': os.path.join(test_dir, "{}".format(genotype)),
        }
    with _open(in_fasta_file) as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            id = seq_record.id
            seq = str(seq_record.seq).upper().translate(trans).replace('-','')
            genotype = sample_mapping[id]
            if id not in dataset[genotype]['test'] and id not in dataset[genotype]['train']:
                continue
            if id in dataset[genotype]['test']:
                dtype = 'test'
            else:
                dtype = 'train'
            fh = open(os.path.join(file_paths[genotype][dtype], "{}.fasta".format(id)),'w')
            fh.write(">{}\n{}\n".format(id,seq))


    return


def run():
    cmd_args = parse_args()
    fasta_file = cmd_args.input_msa
    ref_fasta = cmd_args.ref_fasta
    meta_file = cmd_args.input_meta
    min_genotype_count = cmd_args.min_genotype_count
    max_sample_count = cmd_args.max_samples
    kLen = cmd_args.klen
    delimeter = cmd_args.d
    outdir = cmd_args.outdir
    mask_file = cmd_args.mask
    n_threads = cmd_args.n_threads
    min_length = cmd_args.min_length
    max_ambig = cmd_args.max_ambig
    min_kmer_count = cmd_args.min_kmer_count
    max_homo = cmd_args.max_homo
    stime = time.time()
    if max_homo is None:
        max_homo = int(kLen * 0.2)

    train_proportion = cmd_args.train_proportion

    if not os.path.isdir(outdir):
        os.mkdir(outdir, 0o755)

    analysis_dir = os.path.join(outdir,'_analysis')
    if not os.path.isdir(analysis_dir):
        os.mkdir(analysis_dir, 0o755)

    meta_df = read_metadata(meta_file)
    genotype_count_df = summarize_genotypes(meta_df)
    genotype_count_df['genotype'] = genotype_count_df['genotype'].astype(str)
    genotype_count_df.to_csv(os.path.join(outdir, "genotype_counts.txt"), index=False, sep="\t")
    genotype_count_df = genotype_count_df[genotype_count_df['count'] >= min_genotype_count]
    valid_genotypes = genotype_count_df['genotype'].tolist()
    sample_mapping = read_sample_mapping(meta_df, 'sample_id', 'genotype')
    genotype_counts = {}
    for row in genotype_count_df.itertuples():
        genotype_counts[row.genotype] = row.count

    print("Processing initial alignment")
    stime = time.time()
    seq_data = prep_seqs(fasta_file, sample_mapping, delimeter, analysis_dir, iupac_replacement,max_ambig=max_ambig,min_length=min_length)
    print(time.time() - stime)
    seq_report_file = os.path.join(outdir,"sample.seq.stats.txt")
    sample_report_df = pd.DataFrame.from_dict(seq_data['sample_report'],orient='index')
    sample_report_df.to_csv(seq_report_file,sep="\t",index=False,header=True)
    num_samples = len(seq_data['unique_samples'])

    unique_seq_geno_count = {}
    for md5 in seq_data['unique_samples']:
        sample_id = seq_data['unique_samples'][md5]
        genotype = sample_mapping[sample_id]
        if not genotype in unique_seq_geno_count:
            unique_seq_geno_count[genotype] = 0
        unique_seq_geno_count[genotype]+=1


    print("Calculating genotype consensus seqs")
    stime = time.time()
    genotype_consensus_seqs = {'global':create_consensus_seq(seq_data['global_consensus'])}
    print(time.time() - stime)

    #Create spacer sequence to ensure insertions are preserved
    consensus= seq_data['global_consensus']
    spacer = []
    bases = ['A','T','C','G','N']
    for i in range(0,len(consensus)):
        b = max(consensus[i], key=consensus[i].get)
        if  b == '-':
            for k in bases:
                if consensus[i][k] > 0:
                    b = 'N'
                    break
        spacer.append(b)
    genotype_consensus_seqs['spacer'] = ''.join(spacer)

    tmp_seq_1 = os.path.join(analysis_dir, "seq1.fasta")
    fh = open(tmp_seq_1,'w')
    fh.write(">{}\n{}\n".format('global', genotype_consensus_seqs['global']))
    fh.write(">{}\n{}\n".format('spacer', genotype_consensus_seqs['spacer']))

    for genotype in seq_data['genotype_consensus']:
        if genotype not in valid_genotypes:
            continue
        genotype_consensus_seqs[genotype] = create_consensus_seq(seq_data['genotype_consensus'][genotype])
        fh.write(">{}\n{}\n".format(genotype,genotype_consensus_seqs[genotype]))
    fh.close()
    tmp_out = os.path.join(analysis_dir, "aln.fasta")

    #Align reference to consensus sequences
    print("Aligning reference sequence to genotype consensus alignments")
    stime = time.time()
    mafft_add_seq(ref_fasta, tmp_seq_1, tmp_out, n_threads)
    print(time.time() - stime)
    consensus_seqs = read_fasta(tmp_out)

    ref_seq_dict = read_fasta(ref_fasta)
    ref_id = list(ref_seq_dict.keys())[0]
    ref_aln_lookup = create_aln_pos_from_unalign_pos(consensus_seqs[ref_id])

    #Add ref counts into alignment
    genotype = sample_mapping[ref_id]
    cons = seq_data['genotype_consensus'][genotype]
    seq = consensus_seqs[ref_id]
    for i in range(0,len(seq)):
        cons[i][seq[i]]+=1
    consensus_seqs[genotype] = create_consensus_seq(cons)

    #Get kmer counts for all sequences
    print("Performing kmer counting on all sequences using jellyfish")
    unaligned_file = os.path.join(analysis_dir,"unaligned.fasta")
    stime = time.time()
    if mask_file is not None:
        print("Masking alignment based on supplied file")
        aln_mask = init_mask_regions(len(consensus_seqs['global']), ref_aln_lookup, mask_file)
    else:
        aln_mask = [0] * len(consensus_seqs['global'])
    print(time.time() - stime)
    #Delete gap only columns
    print("Deleting any alignment positions which are only gaps")
    del_mask = init_skip_regions(consensus_seqs)

    #mask positions in each consensus sequence
    consensus_seqs = mask_alignment_dict(consensus_seqs,aln_mask,del_mask,ref_id)

    #Write unaligned consensus sequences for jellyfish counting
    unal_consensus = {}
    consensus_unaln_file = os.path.join(analysis_dir,"consensus.unaln.fasta")
    for seq_id in consensus_seqs:
        unal_consensus[seq_id] = consensus_seqs[seq_id].replace('-','')
    write_fasta_dict(unal_consensus, consensus_unaln_file)

    print("Counting kmers in total dataset")
    stime = time.time()
    jellyfish_result_file = os.path.join(analysis_dir,"jellyfish.txt")
    init_jellyfish_mem = int(num_samples * len(ref_aln_lookup) / 1000000)
    if init_jellyfish_mem == 0:
        init_jellyfish_mem = 1
    jellyfish_mem = "{}M".format(init_jellyfish_mem)
    run_jellyfish_count(unaligned_file, jellyfish_result_file, jellyfish_mem, kLen, n_threads)
    global_kcounts = filter_kmers(get_kmer_counts(jellyfish_result_file), min_kmer_count, num_samples, max_ambig,
                                  max_homo)
    print(time.time() - stime)
    #delete temp file
    os.remove(unaligned_file)
    print("Counting kmers in consensus files")
    stime = time.time()
    jellyfish_cons_file = os.path.join(analysis_dir, "consensus.jellyfish.txt")
    run_jellyfish_count(consensus_unaln_file, jellyfish_cons_file, jellyfish_mem, kLen, n_threads)
    print(time.time() - stime)

    consensus_kcounts = filter_kmers(get_kmer_counts(jellyfish_cons_file), min_kmer_count, num_samples, max_ambig, max_homo)
    needed_kmers = set(list(global_kcounts.keys())) - set(list(consensus_kcounts.keys()))

    num_kmers = len(needed_kmers)
    i = 0
    kmer_index = {}
    for kmer in needed_kmers:
        kmer_index[i] = kmer
        i+=1
    aho = init_automaton_dict(kmer_index)
    print("Performing data reduction for individual genotype sequences")
    stime = time.time()
    #run aho kmer counting on each genotype fasta
    selected_seqs = {}
    num_selected_seqs = 0
    for genotype in valid_genotypes:
        print("{}\t{}".format(genotype,len(needed_kmers)))
        in_fasta_file = os.path.join(analysis_dir,"{}.fasta".format(genotype))
        out_fasta_file = os.path.join(analysis_dir,"{}.selected.fasta".format(genotype))
        result = get_min_seqs_per_genotype(in_fasta_file, out_fasta_file, kmer_index,aho)
        selected_seqs[genotype] =  result['selected_seqs']
        needed_kmers = result['needed_kmers']
        i = 0
        kmer_index = {}
        for kmer in needed_kmers:
            kmer_index[i] = kmer
            i += 1
        aho = init_automaton_dict(kmer_index)
        num_selected_seqs += len(selected_seqs[genotype])
        num_kmers = len(needed_kmers)
    print("Selected {} seqs".format(num_selected_seqs))
    print(time.time() - stime)
    del(result)
    #Write consensus sequences to file
    source_file = os.path.join(analysis_dir,"ref.fasta")
    sourceFH = open(source_file ,'w')
    del(consensus_seqs['spacer'])
    for seq_id in consensus_seqs:
        sourceFH.write(">{}\n{}\n".format(seq_id,consensus_seqs[seq_id]))

    # Merge genotype sequences to file
    target_file = os.path.join(analysis_dir,"alt.fasta")
    mergeFH = open(target_file ,'w')
    for genotype in valid_genotypes:
        in_fasta_file = os.path.join(analysis_dir, "{}.selected.fasta".format(genotype))
        seqs = read_fasta(in_fasta_file)
        if ref_id in seqs:
            del(seqs[ref_id])
        for seq_id in seqs:
            mergeFH.write(">{}\n{}\n".format(seq_id, seqs[seq_id]))
    mergeFH.close()

    print("Aligning selected sequences with reference sequence")
    stime = time.time()
    out_fasta_file = os.path.join(analysis_dir,"aln.fasta")
    mafft_add_seq(source_file,target_file, out_fasta_file, n_threads)
    final_fasta = os.path.join(outdir,"selected_seqs.fasta")
    print(time.time() - stime)
    print("Masking alignment")
    stime = time.time()
    perform_aln_masking_file(out_fasta_file, final_fasta, ref_id, mask_file)
    print(time.time() - stime)

    #create training/testing directory structure
    print("Creating testing/training directory structure")
    train_dir = os.path.join(outdir,'training')
    if not os.path.isdir(train_dir):
        os.mkdir(train_dir, 0o755)

    test_dir = os.path.join(outdir,'test')
    if not os.path.isdir(test_dir):
        os.mkdir(test_dir, 0o755)

    for genotype in selected_seqs:
        genoytpe_dir = os.path.join(train_dir,"{}".format(genotype))
        if not os.path.isdir(genoytpe_dir):
            os.mkdir(genoytpe_dir, 0o755)

        genoytpe_dir = os.path.join(test_dir, "{}".format(genotype))
        if not os.path.isdir(genoytpe_dir):
            os.mkdir(genoytpe_dir, 0o755)

    selected_sample_ids = set()
    datasets = {}
    genotype_train_counts = {}
    for genotype in selected_seqs:
        datasets[genotype] = {'test':set(),'train':set()}
        selected_sample_ids = selected_sample_ids | set(selected_seqs[genotype])
        value = int(train_proportion * genotype_counts[genotype])
        if value > max_sample_count:
            value = max_sample_count
        genotype_train_counts[genotype] = {'needed':value,'selected':0,'test_selected':0}

    pass_seq_df = sample_report_df[sample_report_df['is_len_ok'] == True]
    pass_seq_df = pass_seq_df[pass_seq_df['is_ambig_ok'] == True]
    pass_seq_df = pass_seq_df.sort_values(['len', 'ambig','sample_id'], ascending=[False, True,True])

    print("Dividing dataset into test and train proportions")
    for row in pass_seq_df.itertuples():
        sample_id = row.sample_id
        if sample_id in selected_sample_ids:
            continue
        genotype = row.genotype
        if genotype_train_counts[genotype]['needed'] > genotype_train_counts[genotype]['selected']:
            datasets[genotype]['train'].add(sample_id)
            genotype_train_counts[genotype]['selected']+=1
        elif genotype_train_counts[genotype]['test_selected'] < genotype_train_counts[genotype]['needed'] :
            genotype_train_counts[genotype]['test_selected'] += 1
            datasets[genotype]['test'].add(sample_id)

    print("Writting test and training samples")
    write_individual_files(datasets, sample_mapping, test_dir, train_dir, iupac_replacement, fasta_file)

    #Delete interim files
    #shutil.rmtree(analysis_dir)
    print("{}\t{}".format(num_samples,time.time() - stime))

