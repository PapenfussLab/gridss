[![Release](https://img.shields.io/github/v/release/PapenfussLab/gridss)](https://github.com/PapenfussLab/gridss/releases)
[![Language](http://img.shields.io/badge/language-java-brightgreen.svg)](https://www.java.com/)
[![Language](http://img.shields.io/badge/language-c-brightgreen.svg)]
[![License](https://img.shields.io/badge/license-GPL-blue)]


# VIRUSBreakend - Viral Integration Recognition Using Single Breakends

VIRUSBreakend is a high-speed viral integration detection tool designed to be incorporated in the whole genome sequence piplines with minimal additional cost.
VIRUSBreakend takes on average 1 hour to run on a 100x coverage human sample. Recommended job/VM size is 4 core / 64Gb memory.

This tool is part of GRIDSS - the Genomic Rearrangement IDentification Software Suite.

# Citation

Please cite:

 - VIRUSBreakend: https://doi.org/10.1093/bioinformatics/btab343

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

Run `virusbreakend-build --db virusbreakenddb` to download and generate the reference data.
This download the NCBI taxonomic information, sequences, virushostdb, the `kraken2-build` build process, and generates indexes.
The index is around 54GB in size.
Be aware that the kraken2 build process requires around 150GB of intermediate disk space to download from NCBI and build the index.

`virusbreakend-build` requires:

* An internet connection
* GRIDSS
* samtools
* Kraken2
  * dustmasker (part of the NCBI BLAST+ package)
* [EDirect utilities](https://dataguide.nlm.nih.gov/edirect/install.html)

These can be installed manually, or the BioConda `samtools` `kraken2` and `entrez-direct` packages can be used.

The the generated directory (`virusbreakenddb`) contains the files used in the build process.
The file `virusbreakend.db.virusbreakenddb.tar.gz` contains the minimal subset of the `virusbreakenddb` data required by VIRUSBreakend.


# Running VIRUSBreakend

To run VIRUSbreakend, ensure that the database has been build and run the following command:

```
virusbreakend \
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

## summary.tsv files

|field | meaning|
|----|----|
taxid_genus|NCBI taxonomy ID of genus of viral reference
name_genus|NCBI taxonomy genus name
reads_genus_tree|Number of reads assigned to any virus in that genus by kraken2
taxid_species|NCBI taxonomy ID of species of viral reference
reads_species_tree|Number of reads assigned to that species (and any sub-species) by kraken2
name_species| NCBI taxonomy species name
taxid_assigned|NCBI taxonomy ID of viral reference
name_assigned| NCBI taxonomy name of viral reference
reads_assigned_tree|Number of reads assigned to the viral reference or taxonomic descendents by kraken2
reads_assigned_direct|Number of reads assigned to the viral reference taxid by kraken2
reference|name of viral reference used.
reference_taxid|NCBI taxonomy ID of viral reference used. This can be a child taxid of taxid_assigned.
reference_kmer_count|count of read kmers matching reference used
alternate_kmer_count|count of read kmers matching next best reference candidate
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
QCStatus|QC status of viral integration

Each viral integration should have 2 integration breakpoints (one for the start, one for the end) although there may be more if the integration site is rearranged or less if one side of the  integration site was missed.

## QCStatus

QCStatus can be one of the following:

- ""
The empty string indicates no warning or errors associated with this virus

- LOW_VIRAL_COVERAGE

Less than 10% of the virus has any coverage.
This typically occurs when the viral sequence has a short sequence homology with non-viral sequence.
This record is likely a false positive and should be ignored.

- EXCESSIVE_VIRAL_COVERAGE

Regions of the virus have depth of coverage so high, that GRIDSS was unable to call integrations in these regions.
The exact bounds of excluded regions can be found in the `*.coverage.blacklist.bed` files in the `virusbreakend.working` directory

- ASSEMBLY_DOWNSAMPLED

Viral coverage was sufficiently high that some regions were downsampled in the GRIDSS assembly process.
With the assembly downsampling changes in GRIDSS 2.11.1, viral integration have in these region have probably been correctly called but the QUAL score is lower than it should be.
The exact bounds of excluded regions can be found in the `*.bed.excluded_*.bed` files in the `virusbreakend.working` directory

- CHILD_TAXID_REFERENCE

The viral genome chosen by VIRUSBreakend does not match the NCBI taxonomic classification assigned by Kraken2.
This typically occurs in scenarios such as HPV-45 where Kraken2 assigns reads to the parent taxon due to sequence commonality between HPV-45 and HPV-18.
This is expected behavour for these viral taxa. The taxon assigned `reference_taxid`

- UNCLEAR_TAXID_ASSIGNMENT

Less than 60% of reads assigned to nominal kraken2 assigned taxon were directly assigned (the rest were assigned to child taxa).
The taxonomic assignment may be unclear - possibly due to co-infection by multiple viruses within the genus.

## Ambiguous insertions

The key differentiator of VIRUSBreakend is the ability to detect and classify integration sites in repetative sequences such as centromeres.
Due to the repetative nature of these region, such integration sites cannot be unambigously placed in the host genome.
In such cases, the mapq encoded in the `BEALN` field will be 0 and the field may contain multiple candidicate integration sites.
Integration sites in which the reported position is ambiguous have a `LOW_MAPQ` FILTER applied.

The `INSRM` field contains the repeat sequences identifed in the integration site host sequences.
These annotations can be used to classify ambigous integration sites.

By default, VIRUSBreakend filters out potential integrations with signficant overlap with simple or low complexity repeats since such integrations are typically false positives.
If a more sensitive call set is required, this filter can be removed by rerunning the final `VirusBreakendFilter` steps with the `MINIMUM_REPEAT_OVERLAP` parameter set to any value greater than one (ie require more than 100% repeat overlap).




