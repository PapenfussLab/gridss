# gridss

A high-speed assembly-based next-gen sequencing structural variation caller

# Pre-requisities

To run, gridss the following must be installed:

* Java 1.7 or later
* maven (if building from source)
* NGS aligner. bowtie2 (default), bwa, novoalign, and subread are currently supported

# Building

Pre-compiled binaries are not yet available so gridss must be built from source. Maven is used for build and dependency management which simplifies compile to the following steps:

* git clone https://github.com/d-cameron/idsv
* `mvn package`

If gridss was built successfully, a combined jar containing gridss and all required library located at target/gridss-0.4-SNAPSHOT-jar-with-dependencies.jar

Note: there may or may no still exist some race conditions in the tests causes some tests to fail during the build, but success if run individually, if you encounter this, please raise a github issue. Tests can be skipped by running `mvn package -DskipTests`


# Running

Gridss is built using htsjdk, so is invoked in the same manner as Picard tools utilities. Gridss invokes an external alignment tools at multiple point during processing. When external alignment is required, gridss stops processing and waits for the user to perform the alignment. The gridss processing pipeline performs the following steps:

* gridss: calculate metrics, extract soft clipped reads and discordant read pairs
* aligner: align soft clips
* gridss: generate putative SV assemblies
* aligner: align assemblies
* gridss: call variants

## example/gridss.sh

example/gridss.sh contains an example pipeline of how gridss is invoked.

## libsswjni.so

Due to relatively poor performance of existing Java-based Smith-Waterman alignment packages, gridss incorporates a JNI wrapper to the striped Smith-Waterman alignment library [SSW](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library). lib/libsswjni.so is a precompiled linux version included in gridss If the precompiled version is not compatable with your linux distribution, the wrapper can be compiled by:

* git clone https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library
* make java

If libsswjni.so cannot be loaded, gridss will fall back to a slower java implementation. To ensure libsswjni.so is loaded, add the library using -Djava.library.path, or the LD_LIBRARY_PATH environment variable as per the JNI documentation.

## Parameters

Gridss has a large number of parameters that can be be adjusted. The default parameter set has been tested with paired-end illumina data ranging from 2x36bp to 2x250bp and should give a reasonably good.

Commonly used parameter are listed below, a full list and description of all command-line parameters available with the -h or --help flags.
Debugging/internal usage parameters (for operations such as export of de Bruijn graphs to .gexf) are documented in `Defaults.java`.

### OUTPUT (Required)

Variant calling output file. Can be VCF or BCF.

### REFERENCE (Required)

Reference genome fasta file. Note that gridss caches the entire reference genome in memory so ensure that the working memory available to gridss comfortably exceeds the reference genome size.

### INPUT (Required)

Input libraries. Specify multiple times (ie INPUT=file1.bam INPUT=file2.bam INPUT=file3.bam ) to process multiple libraries together. Gridss considers all reads in each file to come from a single library. BAMs containing read groups from multiple different libraries should be split into per-library

### INPUT_READ_PAIR_MIN_CONCORDANT_FRAGMENT_SIZE, INPUT_READ_PAIR_MAX_CONCORDANT_FRAGMENT_SIZE

Per INPUT overrides of the default concordant fragment size calculation. If your fragment size distribution is unusual, this option can be used to override the.

### INPUT_TUMOUR

Same as INPUT, but used for tumour/normal processing. The somatic P-value (SPV INFO field) is calculated based on total INPUT evidence vs INPUT_TUMOUR evidence.

### INPUT_TUMOUR_READ_PAIR_MIN_CONCORDANT_FRAGMENT_SIZE, INPUT_TUMOUR_READ_PAIR_MAX_CONCORDANT_FRAGMENT_SIZE

Matching INPUT_TUMOUR overrides for INPUT_READ_PAIR_MIN/MAX_CONCORDANT_FRAGMENT_SIZE.

### SCRIPT

When alignment is required to continue processing, gridss will write a script for the alignments required to this file as well as stderr.

### PER_CHR

Flags whether processing should be performed in parallel for each chromosome, or serially for each file.
If your input files are small, or you have a limited number of file handles available, `PER_CHR=false` will
reduce the number of intermediate files created, at the cost of reduced parallelism.

### WORKING_DIR

Directory to write intermediate results directories. By default, intermediate files for each input or output file are written to a subdirectory in the same directory as the relevant input or output file.
If WORKING_DIR is set, all intermediate results are written to subdirectories of the given directory.

### TMP_DIR

This field is a standard Picard tools argument and carries the usual meaning. Temporary files created during processes such as sort are written to this directory.


# Output

Gridss is fundamentally a structural variation breakpoint caller. Variants are output as VCF breakends. Each call is a breakpoint consisting of two breakends, one from location A to location B, and a reciprocal record from location B back to A. Note that although each record fully defines the call, the VCF format required both breakend to be written as separate record.

## Quality score

Gridss calculates quality scores according to the model outlined in [paper].
As gridss does not yet perform multiple test correction or score recalibration, **QUAL scores are vastly overestimated for all variants**.
As a rule of thumb, variants with QUAL >= 1000 and have assembles from both sides of the breakpoint (AS > 0 & RAS > 0) are considered HIGH quality,
variant with QUAL >= 500 but can only be assembled from one breakend (AS > 0 | RAS > 0) are considered MEDIUM quality,
and variants with low QUAL score or lack any supporting assemblies are considered LOW quality.

## Non-standard INFO fields

Gridss writes a number of non-standard VCF fields. These fields are described in the VCF header.


## BEDPE

Gridss supports conversion of VCF to BEDPE format using the VcfBreakendToBedpe utility program included in the gridss jar.

Calling VcfBreakendToBedpe with `INCLUDE_HEADER=true` will include a header containing column names in the bedpe file. These fields match the VCF INFO fields of the same name. For bedpe output, breakend information is not exported and per category totals (such as split read counts) are aggregrated to a single value.








