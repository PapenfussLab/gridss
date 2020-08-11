[![Release](https://img.shields.io/github/v/release/PapenfussLab/gridss)](https://github.com/PapenfussLab/gridss/releases)
[![Build Status](https://travis-ci.org/PapenfussLab/gridss.svg?branch=master)](https://travis-ci.org/PapenfussLab/gridss)
[![Language](http://img.shields.io/badge/language-java-brightgreen.svg)](https://www.java.com/)

# VIRUSBreakend - Viral Integration Recognition Using Single Breakends

- viral integration from BAM file

This tools is part of GRIDSS - the Genomic Rearrangement IDentification Software Suite.

# Citation

Please cite:

 - Preprint for VIRUSBreakend to be uploaded shortly

 - GRIDSS2: harnessing the power of phasing and single breakends in somatic structural variant detection
https://www.biorxiv.org/content/10.1101/2020.07.09.196527v1

 - Kraken2

# Pre-requisites

To run VIRUSBreakend the following must be installed:

* java 1.8 or later
* R 3.6 or later
* samtools
* bwa
* GRIDSS2
* Kraken2

The driver script requires:

* bash
* getopt(1) (part of [util-linux](https://en.wikipedia.org/wiki/Util-linux))

## Reference data setup

- NCBI nodes.dmp https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip

A kraken2 database is required

- Refer to the kraken2 documentation on how to create a kraken2 database.
- The viral library.

Note that building kraken2 requires NCBI `dustbuster`

- *Do not run `kraken2-build --clean`
