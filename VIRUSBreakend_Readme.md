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

 - VIRUSBreakend: https://www.biorxiv.org/content/10.1101/2020.12.09.418731v1

 - GRIDSS2: harnessing the power of phasing and single breakends in somatic structural variant detection
https://www.biorxiv.org/content/10.1101/2020.07.09.196527v1

# Pre-requisites

VIRUSBreakend is part of the GRIDSS software suite. It can be downloaded from https://github.com/PapenfussLab/gridss/releases

All tools used by VIRUSBreakend must be on `PATH` including:
* java
* GRIDSS
* Kraken2
* RepeatMasker
* samtools
* bcftools
* bwa

Set the `GRIDSS_JAR` environment variable to the location of the GRIDSS jar file

## Reference data setup

Run `virusbreakend-build.sh --db virusbreakenddb` to download and generate the reference data.
This download the NCBI taxonomic information, sequences, virushostdb, the `kraken2-build` build process, and generates indexes.
The index is around 7GB in size.
Be aware that the kraken2 build process requires additional memory to build the index.

`virusbreakend-build.sh` requires:

* An internet connection
* GRIDSS
* samtools
* Kraken2
  * dustmasker (part of the NCBI BLAST+ package)

These can be installed manually, or the BioConda samtools and kraken2 packages can be used.

The the generated directory (`virusbreakenddb`) contains the files used in the build process.
The file `virusbreakend.db.virusbreakenddb.tar.gz` contains the minimal subset of the `virusbreakenddb` data required by VIRUSBreakend.


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

VIRUSBreakend outputs:
* A VCF containing the integration breakpoints
* The kraken2 report of the virus(es) for which viral integration was run upon
* Coverage statistics of the virus(es) for which viral integration was run upon

## summary.csv files

|field | meaning|
|----|----|
taxid_genus|NCBI taxonomy ID of genus of viral reference
name_genus|NCBI taxonomy genus name
reads_genus|Number of reads assigned to any virus in that genus by kraken2
taxid_species|NCBI taxonomy ID of species of viral reference
reads_species|Number of reads assigned to that species (and any sub-species) by kraken2
name_species| NCBI taxonomy species name
taxid|NCBI taxonomy ID of viral reference
name| NCBI taxonomy name of viral reference
reads|Number of reads assigned to the viral reference by kraken2
reference|kraken2 name of viral reference contig
#rname|name of adjusted viral reference
startpos|start position of adjusted viral reference. Always 1
endpos|end position of adjusted viral reference. Always equal to the viral contig length
numreads|Number of reads mapped to adjusted viral reference
covbases|Number of bases with at least 1 read mapped
coverage|Percentage of viral sequence with at least one read mapped
meandepth|Mean alignment depth
meanbaseq|Mean base quality of bases mapped to adjusted viral reference
meanmapq|Mean mapping quality of reads mapped to adjusted viral reference
integrations|Number of integration breakpoints found.

Each viral integration should have 2 integration breakpoints (one for the start, one for the end) although there may be more if the integration site is rearranged or less if one side of the  integration site was missed.


## Ambigous insertions

The key differentiator of VIRUSBreakend is the ability to detect and classify integration sites in repetative sequences such as centromeres.
Due to the repetative nature of these region, such integration sites cannot be unambigously placed in the host genome.
In such cases, the mapq encoded in the `BEALN` field will be 0 and the field may contain multiple candidicate integration sites.

The `INSRM` field contains the repeat sequences identifed in the integration site host sequences.
These annotations can be used to classify ambigous integration sites.

By default, VIRUSBreakend filters out potential integrations with signficant overlap with simple or low complexity repeats since such integrations are typically false positives.
If a more sensitive call set is required, this filter can be removed by rerunning the final `VirusBreakendFilter` steps with the `MINIMUM_REPEAT_OVERLAP` parameter set to any value greater than one (ie require more than 100% repeat overlap).




