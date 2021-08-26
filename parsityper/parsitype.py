#!/usr/bin/env python3

import sys
sys.setrecursionlimit(3500)
import multiprocessing as mp

if mp.get_start_method(allow_none=True) != 'spawn':
        mp.set_start_method('spawn', force=True)

tasks = {
    'typer': 'Reconstruct sample genotype(s) from isolate or metagenomic sample',
    'reducer': 'Compress a large MSA into a minimum set of sequences to represent diversity',
    'creator': 'Create a kmer scheme based on labeled data',
    'tuner': 'Train a kmer scheme on labeled genotype data to derive kmer patterns for genotypes',
    'validator': 'Train a kmer scheme on labeled genotype data to derive kmer patterns for genotypes',
    'benchmark': 'Simulate mixed samples and benchmark typing performance',
    'test': 'Test parsityper functionality on a small dataset',
    'version': 'Print version and exit',
}

ordered_tasks = [
    'typer',
    'reducer',
    'creator',
    'tuner',
    'benchmark',
    'validator',
    'test',
    'version'
]


def print_usage_and_exit():
    print('Usage: parsitype <command> [options] <required arguments>', file=sys.stderr)
    print('\nTo get minimal usage for a command use:\nparsitype command', file=sys.stderr)
    print('\nTo get full help for a command use one of:\nparsitype command -h\nparsitype command --help\n', file=sys.stderr)
    print('\nAvailable commands:\n', file=sys.stderr)
    max_task_length = max([len(x) for x in list(tasks.keys())]) + 1
    for task in ordered_tasks:
        print('{{0: <{}}}'.format(max_task_length).format(task), tasks[task], sep=' ', file=sys.stderr)
    sys.exit(0)

def main():

    if len(sys.argv) == 1 or sys.argv[1] in ['-h', '-help', '--help']:
        print_usage_and_exit()

    task = sys.argv.pop(1)

    if task not in tasks:
        print('Task "' + task + '" not recognised. Cannot continue.\n', file=sys.stderr)
        print_usage_and_exit()

    exec('import parsityper.' + task)
    exec('parsityper.' + task + '.run()')

# call main function
if __name__ == '__main__':
    main()