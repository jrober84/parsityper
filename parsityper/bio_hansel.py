from subprocess import Popen, PIPE

class bio_hansel:
    def run_biohansel_directory(scheme_fasta,directory,kmer_file,summary_file,simple_file,min_cov=1,min_frac=0.05,n_threads=1,max_degenerate_kmers=-1):
        if max_degenerate_kmers < 0:
            max_degenerate_kmers = 9999999999999999999
        p = Popen(['hansel',
                   '-s', scheme_fasta,
                   '-D', directory,
                   '-o',summary_file,
                   '-O', kmer_file,
                   '-S', simple_file,
                   '-t',str(n_threads),
                   '--min-kmer-freq',str(min_cov),
                   '--min-kmer-frac',str(min_frac),
                   '--max-degenerate-kmers',str(max_degenerate_kmers),
                   '--force','-vvv'
                   ],
                  stdout=PIPE,
                  stderr=PIPE)
        stdout = str(p.stdout.read())
        stderr = str(p.stderr.read())

        return (stdout, stderr)

    def run_biohansel_single(scheme_fasta,seq_file,kmer_file,summary_file,simple_file,min_cov=8,min_frac=0.05,n_threads=1,max_degenerate_kmers=-1):
        if max_degenerate_kmers < 0:
            max_degenerate_kmers = 9999999999999999999
        p = Popen(['hansel',
                   '-s', scheme_fasta,
                   '-o',summary_file,
                   '-O', kmer_file,
                   '-S', simple_file,
                   '-t',str(n_threads),
                   '--min-kmer-freq',str(min_cov),
                   '--min-kmer-frac',str(min_frac),
                   '--max-degenerate-kmers', str(max_degenerate_kmers),
                   '--force','--verbose'
                   ] + seq_file,
                  stdout=PIPE,
                  stderr=PIPE)
        stdout = str(p.stdout.read())
        stderr = str(p.stderr.read())

        return (stdout, stderr)

    def run_biohansel_pair(scheme_fasta,fwd,rev,kmer_file,summary_file,simple_file,min_cov=8,min_frac=0.05,n_threads=1,max_degenerate_kmers=-1):
        if max_degenerate_kmers < 0:
            max_degenerate_kmers = 9999999999999999999
        p = Popen(['hansel',
                   '-s', scheme_fasta,
                   '-o',summary_file,
                   '-O', kmer_file,
                   '-S', simple_file,
                   '-t',str(n_threads),
                   '--min-kmer-freq',str(min_cov),
                   '--min-kmer-frac',str(min_frac),
                   '--max-degenerate-kmers', str(max_degenerate_kmers),
                   '--force',
                   '-p',fwd,rev
                   ],
                  stdout=PIPE,
                  stderr=PIPE)
        stdout = str(p.stdout.read())
        stderr = str(p.stderr.read())

        return (stdout, stderr)
