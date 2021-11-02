from parsityper.ext_tools import run_command

def run_canu_correct(reads,prefix,out_dir,genome_size,min_length=500,minOverlapLength=500,corOutCoverage=1000,n_threads=1):
    cmd_args = {'-t ':n_threads,
                '-p ':prefix,
                'stopOnLowCoverage=':'False',
                'genomeSize=':genome_size,
                'minReadLength=':min_length,
                'minOverlapLength':minOverlapLength,
                'corOutCoverage':corOutCoverage,
                '-d ':out_dir,
                '-nanopore ':reads}
    cmd = "canu -correct useGrid=False {}".format((" ".join(f'{k}{v}' for k,v in cmd_args.items())))
    return run_command(cmd)