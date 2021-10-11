# GRIDSS Quick start guide

What do you want to do?

- [Install GRIDSS](#install-gridss)
  * [From conda](#from-conda)
  * [From github releases](#from-github-releases)
- [Call structural variants](#call-structural-variants)
  * [Which blacklist file should I use?](#which-blacklist-file-should-i-use-)
- [Call somatic structural variants](#call-somatic-structural-variants)
- [Find viral integration sites in a human genome](#find-viral-integration-sites-in-a-human-genome)
- [Find viral integration sites in a non-human genome](#find-viral-integration-sites-in-a-non-human-genome)
- [Re-call structural variants from another SV caller](#re-call-structural-variants-from-another-sv-caller)
- [Call structural variants in a small region of the genome:](#call-structural-variants-in-a-small-region-of-the-genome-)
- [Annotate single breakend variant calls](#annotate-single-breakend-variant-calls)
  * [RepeatMasker](#repeatmasker)
  * [Kraken2](#kraken2)
- [Optimise GRIDSS execution](#optimise-gridss-execution)
  * [setupreference](#setupreference)
  * [preprocess](#preprocess)
  * [assembly](#assembly)
  * [Variant Calling](#variant-calling)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>

# Install GRIDSS

## From conda

- Install [BioConda](https://bioconda.github.io/)

- Create and use a GRIDSS environment

```
conda create -n gridss gridss
conda activate gridss
```

## From github releases

Download the latest release from https://github.com/PapenfussLab/gridss/releases.

Install all [GRIDSS dependencies](https://github.com/PapenfussLab/gridss#pre-requisites)


# Call structural variants

This can be done in a single command:
```
gridss \
  -r reference.fa \
  -j gridss-2.X.Y-gridss-jar-with-dependencies.jar \
  -o output.vcf \
  -b blacklist.bed \
  input1.bam \
  input2.bam
```


If your input files are aligned with `bwa mem` or another aligner that reports split read alignments using the `SA` tag, then runtime can be reduced by specifying `--skipsoftcliprealignment`.

## Which blacklist file should I use?

We recommend using the ENCODE Blacklist files when using human data:

|Genome|Blacklist|
|-|-|
|hg19 (NCBI chr notation)|ENCFF001TDO.bed|
|hg38 (UCSC chr notation)|ENCFF356LFX.bed|

A mirror of these files can be found [here](https://github.com/PapenfussLab/gridss/tree/master/example)

# Call somatic structural variants

First perform joint SV calling on the tumour and normal.
Make sure that you specify the normal bam then the tumour bam.

```
gridss \
  -r reference.fa \
  -j gridss-2.X.Y-gridss-jar-with-dependencies.jar \
  -o all_calls.vcf \
  -b blacklist.bed \
  normal.bam \
  tumour.bam
```

Then filter the output VCF down to just the somatic SVs:

```
gridss_somatic_filter \
  --pondir refdata/hg19/dbs/gridss/pon3792v1/ \
  --input all_calls.vcf \
  --output high_confidence_somatic.vcf.gz \
  --fulloutput high_and_low_confidence_somatic.vcf.gz \
  --scriptdir $(dirname $(which gridss_somatic_filter)) \
  -n 1 \
  -t 2
```

If you are using hg19 or hg38, then a PON based on Dutch samples is available from https://resources.hartwigmedicalfoundation.nl/ in the HMFTools-Resources/GRIDSS directory.

Note that for the PON filtering to work:

- Both the `.bedpe` and `.bed` must be present in the `--pondir` directory.
- They must be named `gridss_pon_breakpoint.bedpe` and `gridss_pon_single_breakend.bed` respectively.
- The chromosome naming conventions (NCBI-style `1` vs UCSC-style `chr1`) must match the convention used by your reference genome.

# Find viral integration sites in a human genome

Download and extract the pre-built VIRUSBreakend database

```
wget https://virusbreakend.s3.us-east-2.amazonaws.com/virusbreakenddb_20210401.tar.gz
tar zxvf https://virusbreakend.s3.us-east-2.amazonaws.com/virusbreakenddb_20210401.tar.gz
```

Then run VIRUSBreakend:

```
virusbreakend \
  -r host_reference.fa \
  -j gridss-2.X.Y-gridss-jar-with-dependencies.jar \
  -o output.vcf \
  --db virusbreakenddb_20210401 \
  host_aligned_input.bam
```

# Find viral integration sites in a non-human genome

Unfortunately, the prebuild Kraken2 database won't work for non-human genomes as human is the only host genome supplied by Kraken2.
You'll need to build your own.

First build the human Kraken2 reference:

```
virusbreakend-build --db my_custom_virusbreakend_database
```

As VIRUSBreakend does not officially support non-human hosts, you'll need to tell VIRUSBreakend to treat your host as human.

Edit your host `.fa` file and prepend `kraken:taxid|9606|` to every contig name.
This means that `>chr1` becomes `>kraken:taxid|9606|chr1`, and so on.
```
sed 's/^>/>kraken:taxid|9606|/' < my_host_reference.fa > my_host_kraken9606_reference.fa
```

Once you've edited your reference, add it to the VIRUSBreakend database, and rebuild.

```
kraken2-build --add-to-library my_host_kraken9606_reference.fa --db my_custom_virusbreakend_database
kraken2-build $kraken2buildargs --build --db my_custom_virusbreakend_database
```

With the custom database built, you can now run `virusbreakend` against your custom database:

```
virusbreakend \
  -r my_host_reference.fa \
  -j gridss-2.X.Y-gridss-jar-with-dependencies.jar \
  -o output.vcf \
  --db my_custom_virusbreakend_database \
  host_aligned_input.bam
```


# Re-call structural variants from another SV caller

First extract the fragments overlapping the region of interest:

```
gridss_extract_overlapping_fragments \
  --targetvcf othercaller.vcf \
  -j gridss-2.X.Y-gridss-jar-with-dependencies.jar \
  input.bam
```

Then run GRIDSS on the tageted bam:

```
gridss \
  -r reference.fa \
  -j gridss-2.X.Y-gridss-jar-with-dependencies.jar \
  -o output.vcf \
  input.bam.targeted.bam
```

# Call structural variants in a small region of the genome:

First extract the fragments overlapping the region of interest:

```
gridss_extract_overlapping_fragments \
  --targetbed region_of_interest.bed \
  -j gridss-2.X.Y-gridss-jar-with-dependencies.jar \
  input.bam
```

Then run GRIDSS on the tageted bam:

```
gridss \
  -r reference.fa \
  -j gridss-2.X.Y-gridss-jar-with-dependencies.jar \
  -o output.vcf \
  input.bam.targeted.bam
```

# Annotate single breakend variant calls

## RepeatMasker

To annotate single breakends (and sequences inserted in breakpoints) with RepeatMasker annotations, ensure RepeatMasker is on `PATH` and run:

```
gridss_annotate_vcf_repeatmasker \
  -j gridss-2.X.Y-gridss-jar-with-dependencies.jar \
  -o annotated.vcf \
  input.vcf
```

## Kraken2

Kraken2 annotation can be useful to determine if single breakend variant calls are the result of viral integration.
To annotate single breakends (and sequences inserted in breakpoints) with a species classification, ensure Kraken2 is on `PATH` and run:

```
gridss_annotate_vcf_repeatmasker \
  -j gridss-2.X.Y-gridss-jar-with-dependencies.jar \
  -o annotated.vcf \
  --kraken2db standard
  input.vcf
```

Pre-build Kraken2 database can be downloaded from the [Kraken2 AWS Public Dataset Program](https://benlangmead.github.io/aws-indexes/k2).
Alternatively, the VIRUSBReakend database can also be used.


# Optimise GRIDSS execution

If running GRIDSS in a cluster compute environment, additional parallelisation options are available.

Note that `--workingdir` must be the same directory for all steps.

## setupreference

Firstly, a once-off setup step for the reference genome is required.
This only needs to be run once per reference genome.

```
gridss \
  -r reference.fa \
  -j gridss-2.X.Y-gridss-jar-with-dependencies.jar \
  -s setupreference
```

## preprocess

Pre-processing is done per input file and has good scaling to 4 cores.

Peak memory usage occurs either during bwa alignment (5Gb + bwa index size) `--skipsoftcliprealignment` is not specified, or during sort (5Gb + 768Mb per thread).

```
gridss \
  -r reference.fa \
  -j gridss-2.X.Y-gridss-jar-with-dependencies.jar \
  -s preprocess \
  -t 4 \
  input.bam
```

## assembly

Assembly is performed in parallel on 10Mbase chunks (chunk size can be adjusted by setting a different value through a custom `gridss.properties` and the `--configuration` parameter).

Assembly scales to 8 cores with peak memory usage of 31Gb and can be distributed across multiple nodes.

Running multiple assembly jobs requires specifying `--jobnodes` and a unique 0-based `--jobindex` for each job.

For example, to distribute assembly of paired tumour/normal sequencing into 3 jobs, the jobs would be:

```
gridss -r reference.fa -j gridss-2.X.Y-gridss-jar-with-dependencies.jar \
  -t 8 \
  -s assemble \
  -a assembly.bam \
  --jobnodes 3 \
  --jobindex 0 \
  normal.bam tumour.bam
```
```
gridss -r reference.fa -j gridss-2.X.Y-gridss-jar-with-dependencies.jar \
  -t 8 \
  -s assemble \
  -a assembly.bam \
  --jobnodes 3 \
  --jobindex 1 \
  normal.bam tumour.bam
```
```
gridss -r reference.fa -j gridss-2.X.Y-gridss-jar-with-dependencies.jar \
  -t 8 \
  -s assemble \
  -a assembly.bam \
  --jobnodes 3 \
  --jobindex 2 \
  normal.bam tumour.bam
```

A final assembly job is required to merge the per-chunk outputs together.
This job fast, I/O-bound, and is best combined with the variant calling job:

## Variant Calling

The final variant calling step also uses 8 cores and 31Gb:

```
gridss -r reference.fa -j gridss-2.X.Y-gridss-jar-with-dependencies.jar \
  -t 8 \
  -a assembly.bam \
  -s assemble,call \
  -o output.vcf \
  normal.bam tumour.bam
```


