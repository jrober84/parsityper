|logo|

|conda| |nbsp| |pypi| |nbsp|  |rtd| |nbsp| |license|


======  ===========
Master  |ci-master|
Dev     |ci-dev|
======  ===========

.. |logo| image:: logo.png
    :target: https://github.com/jrober84/parsityper


Introduction
============
**ALPHA version** Not all biological sample types are true isolates whether that is due to the uncultivability of the organism or the analysisof complex biological samples such as food matrices or primary clinical samples. Occam's razor or law of parsimony — "entities should not be multiplied without necessity". Parsityper applies k-mer genotyping schemes to assembled or raw sequence data to find the smallest number of compatible genotypes which explain the provided data. In cases where a k-mer target cannot be assigned to a genotype it will be labelled as “unknown”.


Citation
========

Requirements and Dependencies
=============================
- bio_hansel >= 2.6.1
- Python_ (>=v3.6)
    - numpy_ >=1.12.1
    - pandas_ >=1.1.3
    - pyahocorasick_ >=1.1.6
    - attrs_
    - biopython >= 1.7.0
    - scipy >= 1.52



Installation
============

With pip_ from Github
---------------------

Install the latest master branch version directly from Github:

.. code-block:: bash

    pip install git+https://github.com/jrober84/parsityper.git

Usage
=====
If you run ``parsityper``, you should see the following usage statement:

.. code-block::

    Usage: parsitype <command> [options] <required arguments>

    To get minimal usage for a command use:
    parsitype command

    To get full help for a command use one of:
    parsitype command -h
    parsitype command --help


    Available commands:

    typer    Reconstruct sample genotype(s) from isolate or metagenomic sample [alpha]
    creator  Create a kmer scheme based on labeled data [ alpha]
    trainer  Train a kmer scheme on labeled genotype data to derive kmer patterns for genotypes [not implemented]
    test     Test parsityper functionality on a small dataset [ not implemented]
    version  Print version and exit

**Typer**
typer --mode batch --type multi --scheme SARS-COV-2_v1 --data_dir {input_directory} --outdir {result_directory}

Optional paramters:
prefix - change the prefixe of result files
mode - 'batch' or 'single' [**not implemented**] But will trigger run level analysis between samples to look for issues
type - 'multi' or 'single' Toggles QA/QC for metagenomic or isolate samples respectively, should set to what is expectation regarding the samples
primers - default (none) There are several inbuilt primer schemes 'arctic_v1', 'arctic_v2', 'arctic_v3', 'freed_v1', 'resende_v1' which can be used to subtract heterozygous sites which overlap primers from quality analysis
min_cov - default (50X) Not recommended to change right now as this seems to work well for sars-cov-2 datasets
min_cov_frac - default (0.05) Illumina data this is the lower limit, Nanopore data (0.2) should be the lower limit. Increase this to reduce potential false positive minor variants
max_mixed_sites - default(10) Illumina data, most datasets have lower than this. Nanopore data set to 15 to flag mixed samples
max_missing_sites - default(1400) Flags samples with poor k-mer target coverage. This value works for both Illumina and Nanopore
n_threads - very minor performance difference on practical sized dataset with more than 32 threads


Output files:
{prefix}.sample_composition.report.summary.txt - Sample summary report of mutations and genotypes present in the sample
{prefix}.sample_composition.report.detailed.txt - Mutation associations to genotype for each sample
{prefix}.run_kmer_dendrogram.png - If more than one sample is analysed then the jaccard similarity is estimated between samples to produce a dendrogram
{prefix}.primers.report.txt - Detected kmers in primers (sample_id,target,seq,is_positive,freq)
{prefix}.positive_control.report.txt - Detected sample kmers in supplied positive control (sample_id,target,seq,is_positive,freq)
{prefix}.negative_control.report.txt - Detected sample kmers in supplied not template control (sample_id,target,seq,is_positive,freq)
{prefix}.run_composition.report.txt - K-mer targets and the number of samples with heterozygous calls



Legal
=====

Copyright Government of Canada 2021

Written by: National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.


Contact
=======

**James Robertson**: james.robertson@canada.ca