import shutil

from parsityper.ext_tools import run_command
import re, os

def kmc_command(fastq_read,out_file, tmp_file,mode='fastq',freq=10,kmer_length=21,n_threads=1):
    cmd_args = {
        '-sm': '',
        '-m':12,
        '-k': kmer_length,
        '-t': n_threads,
    }
    if mode == 'fastq':
        cmd_args['-fq'] = ''
        cmd_args['-ci'] = freq
    else:
        cmd_args['-fm']= ''
        cmd_args['-ci'] = 1

    cmd = "kmc {} {} {} {}".format((" ".join(f'{k}{v}' for k,v in cmd_args.items())),fastq_read,out_file,tmp_file)
    return run_command(cmd)


def kmc_summary(fastq_read,out_file, tmp_dir, mode,freq=10,kmer_length=21,n_threads=1):
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir, 0o755)

    data = {'est_genome_size':0,'num_unique_kmers':0,'num_counted_unique_kmers':0,'num_reads':0}
    (stdout,stderr) = kmc_command(fastq_read,out_file, tmp_dir,mode,freq,kmer_length,n_threads)
    m = re.findall("No. of unique counted k-mers.+\s+(\d+)",stdout)
    if len(m) != 0:
        data['est_genome_size'] = int(m[0])
        data['num_counted_unique_kmers'] = data['est_genome_size']
    m = re.findall("No. of unique k-mers.+\s+(\d+)",stdout)
    if len(m) != 0:
        data['num_unique_kmers'] = int(m[0])
    m = re.findall("Total no. of reads.+\s+(\d+)",stdout)
    if len(m) != 0:
        data['num_reads'] = int(m[0])

    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    if os.path.isfile("{}.kmc_pre".format(out_file)):
        os.remove("{}.kmc_pre".format(out_file))
        os.remove("{}.kmc_suf".format(out_file))

    return data