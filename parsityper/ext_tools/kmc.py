from parsityper.ext_tools import run_command
import re

def kmc_command(fastq_read,out_file, tmp_file,freq=10,kmer_length=21,n_threads=1):
    cmd_args = {
        '-sm': '',
        '-m': freq,
        '-k': kmer_length,
        '-t': n_threads,
    }
    cmd = "kmc {} {} {}".format((" ".join(f'{k}{v}' for k,v in cmd_args.items())),fastq_read,out_file,tmp_file)
    return run_command(cmd)


def kmc_summary(fastq_read,out_file, tmp_file,freq=10,kmer_length=21,n_threads=1):
    data = {'est_genome_size':0,'num_unique_kmers':0,'num_counted_unique_kmers':0,'num_reads':0}
    (stdout,stderr) = kmc_command(fastq_read,out_file, tmp_file,freq,kmer_length,n_threads)
    m = re.findall("No. of unique counted k-mers.+\s+(\d+)",stdout)
    if len(m) == 0:
        data['est_genome_size'] = int(m[0])
        data['num_counted_unique_kmers'] = data['est_genome_size']
    m = re.findall("No. of unique k-mers.+\s+(\d+)",stdout)
    if len(m) == 0:
        data['num_unique_kmers'] = int(m[0])
    m = re.findall("Total no. of reads.+\s+(\d+)",stdout)
    if len(m) == 0:
        data['num_reads'] = int(m[0])
    return data