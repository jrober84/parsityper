import os.path
import json
from parsityper.ext_tools import run_command



def run_fastp(read_set,out_dir,out_prefix,min_read_len=0,max_read_len=0,trim_front_bp=0,trim_tail_bp=0,report_only=True,dedup=False,merge_reads=False,n_threads=1):
    json = os.path.join(out_dir,"{}.json".format(out_prefix))
    html = os.path.join(out_dir, "{}.html".format(out_prefix))
    out1 = os.path.join(out_dir, "{}_1.fastq".format(out_prefix))
    out2 = os.path.join(out_dir, "{}_2.fastq".format(out_prefix))
    merged_out = os.path.join(out_dir, "{}.merged.fastq".format(out_prefix))
    cmd_args = {'-j ':json, '-h ':html, '-w ':n_threads}
    cmd_args['-i '] = read_set[0]
    cmd_args['-f '] = trim_front_bp
    cmd_args['-t '] = trim_tail_bp
    cmd_args['-l '] = min_read_len
    cmd_args['--length_limit '] = max_read_len
    if dedup:
        cmd_args['-D '] = ''
    if not report_only:
        cmd_args['-o '] = out1
    if len(read_set) == 2:
        cmd_args['-I '] = read_set[1]
        cmd_args['-F '] = trim_front_bp
        cmd_args['-T '] = trim_tail_bp
        if not report_only:
            cmd_args['-O '] = out2
        if merge_reads:
            cmd_args['-m'] = ''
            cmd_args['--merged_out '] = merged_out



    cmd = "fastp {}".format((" ".join(f'{k}{v}' for k,v in cmd_args.items())))
    (stdout,stderr) = run_command(cmd)
    return process_json(json)

def process_json(json_path):
    with open(json_path) as json_file:
        return json.load(json_file)

    return {}
