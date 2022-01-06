Welcome to Parsityper Read the Docs!
=====================================

Parsityper (parsimonious genotyper) is designed as a suite of tools designed to build, validate and utilize k-mer genotyping schemes based on predetermined genotypes with an emphasis on public health diagnostic laboratories.  This includes an emphasis on “batch” analysis with positive and negative control samples to provide insights into potential contamination or technical issues with the sequence data in addition to standard individual sequence QC metrics. The driving principal behind the typing module is to find the most parsimonious (simplest) genotype composition for the data using k-mers which are informative for each genotype. The typing module is compatible consensus fasta or raw fastq data from any sequencing technology (Illumina, Ion Torrent, Nanopore, PacBio) and is designed to work with both isolate and complex sample types (i.e. waste-water metagenomes) using standard WGS or amplicon-based approaches. For amplicon-based sequences, there is support for removing primer sequences from reads prior to genotype analysis.

Code is available on GitHub under https://github.com/jrober84/parsityper.

* :ref:`user-docs`
* :ref:`legal`
* :ref:`contact`


.. _user-docs:

.. toctree::
   :maxdepth: 2
   :caption: User Documentation
   user-docs/installation
   user-docs/usage-composer
   user-docs/usage-creator
   user-docs/usage-tuner
   user-docs/usage-validator
   user-docs/usage-benchmark
   user-docs/usage-typer
   user-docs/Tutorial
   user-docs/genotyping_schemes
   user-docs/degenerate_base_expansion
   user-docs/output
   user-docs/parameters
   user-docs/command-line

.. _legal:

.. toctree::
   :maxdepth: 2
   :caption: Legal

   legal

.. _contact:

.. toctree::
   :maxdepth: 2
   :caption: Contact

   contact

.. toctree::
   :maxdepth: 2
   :caption: poop