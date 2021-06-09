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
Not all biological sample types are true isolates whether that is due to the unculurability of the organism or the analysis
of complex biological samples such as food matricies or primary clinical samples.
Occam's razor or law of parsimony — "entities should not be multiplied without necessity"



Citation
========

Requirements and Dependencies
=============================

Installation
============

Usage
=====
If you run ``parsityper``, you should see the following usage statement:

.. code-block::

    usage: hansel [-h] [-s SCHEME] [--scheme-name SCHEME_NAME]
                  [-M SCHEME_METADATA] [-p forward_reads reverse_reads]
                  [-i fasta_path genome_name] [-D INPUT_DIRECTORY]
                  [-o OUTPUT_SUMMARY] [-O OUTPUT_KMER_RESULTS]
                  [-S OUTPUT_SIMPLE_SUMMARY] [--force] [--json]
                  [--min-kmer-freq MIN_KMER_FREQ] [--min-kmer-frac MIN_KMER_FRAC]
                  [--max-kmer-freq MAX_KMER_FREQ]
                  [--low-cov-depth-freq LOW_COV_DEPTH_FREQ]
                  [--max-missing-kmers MAX_MISSING_KMERS]
                  [--min-ambiguous-kmers MIN_AMBIGUOUS_KMERS]
                  [--low-cov-warning LOW_COV_WARNING]
                  [--max-intermediate-kmers MAX_INTERMEDIATE_KMERS]
                  [--max-degenerate-kmers MAX_DEGENERATE_KMERS] [-t THREADS] [-v]
                  [-V]
                  [F [F ...]]

    BioHansel version 2.5.1: Subtype microbial genomes using SNV targeting k-mer subtyping schemes.

    Built-in schemes:

    * heidelberg:  Salmonella enterica spp. enterica serovar Heidelberg
    * enteritidis: Salmonella enterica spp. enterica serovar Enteritidis
    * typhimurium: Salmonella enterica spp. enterica serovar Typhimurium
    * typhi:       Salmonella enterica spp. enterica serovar Typhi
    * tb_lineage:  Mycobacterium tuberculosis

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