[![Release](https://img.shields.io/github/v/release/PapenfussLab/gridss)](https://github.com/PapenfussLab/gridss/releases)
[![Build Status](https://travis-ci.org/PapenfussLab/gridss.svg?branch=master)](https://travis-ci.org/PapenfussLab/gridss)
[![Language](http://img.shields.io/badge/language-java-brightgreen.svg)](https://www.java.com/)
[![Language](http://img.shields.io/badge/language-c-brightgreen.svg)]
[![License](https://img.shields.io/badge/license-GPL-blue)]


# VIRUSBreakend - Viral Integration Recognition Using Single Breakends

VIRUSBreakend is a high-speed viral integration detection tool designed to be incorporated in the whole genome sequence piplines with minimal additional cost.
VIRUSBreakend takes on average 1 hour to run on a 100x coverage human sample (4 core c2-standard-4 google compute instance) at a cost of USD$0.21 per sample (USD$0.06 using preemptible instances).

This tool is part of GRIDSS - the Genomic Rearrangement IDentification Software Suite.

# Citation

Please cite:

 - Preprint for VIRUSBreakend to be uploaded shortly

 - GRIDSS2: harnessing the power of phasing and single breakends in somatic structural variant detection
https://www.biorxiv.org/content/10.1101/2020.07.09.196527v1

# Pre-requisites

To run VIRUSBreakend the following must be installed:

* java 1.8 or later
* R 3.6 or later
* samtools
* bwa
* Kraken2
* RepeatMasker
* htslib 1.10
* GRIDSS

The driver script requires:

* bash
* getopt(1) (part of [util-linux](https://en.wikipedia.org/wiki/Util-linux))

Once 
* Ensure GRIDSS, Kraken2, RepeatMasker, samtools and bwa are on `PATH`
* Set the `GRIDSS_JAR` environment variable to the location of the GRIDSS jar file


## gridsstools

Performance-critical steps in VIRUSBreakend are implemented in C using htslib.
A precompiled version of `gridsstools` is included as part of GRIDSS releases.
If this precompiled version does not run on your system you will need to build it from source.

To build `gridsstools` from source, download and install htslib then run:

```
git clone http://github.com/PapenfussLab/gridss/
cd gridss/src/main/c/gridsstools/
autoheader
autoconf
./configure && make all
```

## Reference data setup

Run `virusbreakend-build.sh virusbreakenddb` to download and generate the reference data.

`virusbreakend-build.sh` requires:

* An internet connection
* samtools
* Kraken2
  * dustmasker (part of the NCBI BLAST+ package)

These can be installed manually, or the BioConda samtools and kraken2 packages can be used.

The the generated file `virusbreakend.db.virusbreakenddb.tar.gz` contains the viral reference data used by VIRUSBreakend.

# Running VIRUSBreakend

To run VIRUSbreakend, ensure that the database has been build and run the following command:

```
virusbreakend.sh \
	--kraken2db virusbreakenddb \
	--output sample.virusbreakend.vcf \
	--reference host_reference.fa \
	sample.bam
```

# Output

The output format is a VCF file containing the location of single breakend from the viral sequence.
The integration location in the host is encoded in the `BEALN` field.
Note that depending on the host alignment and single breakend orientations, the integration position will be at either the start or end of the `BEALN` alignment position.

In future versions, this is likely to be replaced by a more readable breakpoint `BND` notation.

## Ambigous insertions

The key differentiator of VIRUSBreakend is the ability to detect and classify integration sites in repetative sequences such as centromeres.
Due to the repetative nature of these region, such integration sites cannot be unambigously placed in the host genome.
In such cases, the mapq encoded in the `BEALN` field will be 0 and the field may contain multiple candidicate integration sites.

The `INSRM` field contains the repeat sequences identifed in the integration site host sequences.
These annotations can be used to classify ambigous integration sites.

By default, VIRUSBreakend filters out potential integrations with signficant overlap with simple or low complexity repeats since such integrations are typically false positives.
If a more sensitive call set is required, this filter can be removed by rerunning the final `VirusBreakendFilter` steps with the `MINIMUM_REPEAT_OVERLAP` parameter set to any value greater than one (ie require more than 100% repeat overlap).




