import os
import pandas as pd
from parsityper.ext_tools import run_command

def mash_sketch_cmd(input_files,out_prefix,is_Reads,min_copies=1,sketch_size=1000,individual_sketch=False,kmer_length=21,seed=42,n_threads=1):
    tmp_file_path = os.path.join(out_prefix,".tmpfile")
    tmp_file = open(tmp_file_path,'w')
    for f in input_files:
        tmp_file.write("{}\n".format(f))
    tmp_file.close()
    cmd_args = {
        '-l':tmp_file_path,
        '-o': out_prefix,
        '-s':sketch_size,
        '-k':kmer_length,
        '-S':seed,
        '-p':n_threads,
    }
    if individual_sketch:
        cmd_args['-i'] = ''
    if is_Reads:
        cmd_args['-r'] = ''
        cmd_args['-m'] = min_copies

    cmd = "mash screen {}".format((" ".join(f'{k} {v}' for k,v in cmd_args.items())))
    return run_command(cmd)

def mash_dist_cmd(reference_sketch,query_sketch,max_pvalue=1,max_dist=1,sketch_size=1000,individual_sketch=False,table_fmt=False,kmer_length=21,seed=42,n_threads=1):
    cmd_args = {
        '-s':sketch_size,
        '-k':kmer_length,
        '-S':seed,
        '-p':n_threads,
        '-v':max_pvalue,
        '-d':max_dist
    }
    if individual_sketch:
        cmd_args['-i'] = ''
    if table_fmt:
        cmd_args['-t'] = ''

    cmd = "mash dist {} {} {}".format((" ".join(f'{k} {v}' for k,v in cmd_args.items())),reference_sketch,query_sketch)
    return run_command(cmd)

def mash_screen_cmd(reference_sketch,query_file,out_file,max_pvalue=1,max_dist=1,winner_take_all=True,n_threads=1):
    cmd_args = {
        '-p':n_threads,
        '-v':max_pvalue,
        '-d':max_dist
    }
    if winner_take_all:
        cmd_args['-w'] = ''

    cmd = "mash screen {} {} {}".format((" ".join(f'{k} {v}' for k,v in cmd_args.items())),query_file,reference_sketch)
    return run_command(cmd)

def mash_paste_cmd(out_prefix,paths_file):
    cmd = "mash paste -l {} {}".format(out_prefix,paths_file)
    return run_command(cmd)

def mash_sample_comparisons(sampleManifest,mash_dir,min_freq=10,k=21,s=1000,seed=42,n_threads=1):
    #Sketch files
    sketches = {}
    for sampleID in sampleManifest:
        if len(sampleManifest[sampleID]['processed_reads'])> 0:
            read_set = sampleManifest[sampleID]['processed_reads']
        else:
            read_set = sampleManifest[sampleID]['raw_seq_files']
        sketch_path = os.path.join(mash_dir,sampleID)
        is_Reads = False
        min_copies = 1
        if sampleManifest[sampleID]['seq_type'] == 'fastq':
            is_Reads = True
            min_copies = min_freq
        sketches[sampleID] = "{}.msh".format(sketch_path)
        (stdout,stderr) = mash_sketch_cmd(read_set, sketch_path, is_Reads, min_copies=min_copies, sketch_size=s, individual_sketch=False,
                        kmer_length=k, seed=seed, n_threads=n_threads)

    #combine sketches
    sketch_manifest = os.path.join(mash_dir,"sketches.txt")
    sketch_prefix = os.path.join(mash_dir,"merged.sketches")
    combined_sketch_file = "{}.msh".format(sketch_prefix)
    fh  = open(sketch_prefix,'w')
    for sampleID in sketches:
        fh.write("{}\n".format(sketches[sampleID]))
    fh.close()
    (stdout, stderr) = mash_paste_cmd(sketch_prefix,sketch_manifest)

    #calculate distances
    mash_dist_mat = os.path.join(mash_dir,"distance.matrix.txt")
    (stdout, stderr) = mash_dist_cmd(combined_sketch_file,
                                     combined_sketch_file,
                                     max_pvalue=1,
                                     max_dist=1,
                                     sketch_size=s,
                                     individual_sketch=False,
                                     table_fmt=True,
                                     kmer_length=k,
                                     seed=seed,
                                     n_threads=n_threads)
    fh = open(mash_dist_mat, 'w')
    fh.write(mash_dist_mat)
    fh.close()

    df = pd.read_csv(mash_dist_mat,sep="\t",header=0)
    sampleIDs = df['#query'].tolist()
    df = df.drop(['#query'])
    return {'samples':sampleIDs,'matrix':df.to_numpy()}
