from parsityper.ext_tools import run_command

def run_lighter(read_set,out_dir,genome_size_alpha,genome_size,n_threads=1):
    cmd_args = {'-r':read_set[0], '-r':read_set[0], '-t':n_threads,'-k':genome_size_alpha,'-K':genome_size,'-od':out_dir}
    cmd = "lighter {}".format((" ".join(f'{k}{v}' for k,v in cmd_args.items())))
    return run_command(cmd)

