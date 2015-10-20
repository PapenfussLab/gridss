# gridss

A high-speed next-gen sequencing structural variation caller.
gridss calls variants based on alignment-guided positional de Bruijn graph breakpoint assembly, split read, and read pair evidence.

# Pre-requisities

To run, gridss the following must be installed:

* Java 1.8 or later
* maven (if building from source)
* NGS aligner. bowtie2 (default), bwa, novoalign, and subread are currently supported

## Building from source

Maven is used for build and dependency management which simplifies compile to the following steps:

* git clone https://github.com/PapenfussLab/gridss
* `mvn package -DskipTests`

If gridss was built successfully, a combined jar containing gridss and all required library located at target/gridss-*-jar-with-dependencies.jar

# Running

Pre-compiled binaries are available at https://github.com/PapenfussLab/gridss/releases.

Gridss is built using htsjdk, so is invoked in the same manner as Picard tools utilities. Gridss invokes an external alignment tools at multiple point during processing. When external alignment is required, gridss stops processing and waits for the user to perform the alignment.
The gridss processing pipeline performs the following steps:

* gridss: calculate metrics, extract soft clipped reads and discordant read pairs
* aligner: align soft clips
* gridss: generate putative SV assemblies
* aligner: align assemblies (x3 to handle compound breakpoints)
* gridss: call variants

## example/gridss.sh

example/gridss.sh contains an example pipeline of how gridss is invoked.

## Parameters

Gridss has a large number of parameters that can be be adjusted. The default parameter set has been tested with paired-end illumina data ranging from 2x36bp to 2x250bp and should give a reasonably good.

Commonly used parameter are listed below, a full list and description of all command-line parameters available with the -h or --help flags.
Debugging/internal usage parameters (for operations such as export of de Bruijn graphs to .gexf) are documented in `Defaults.java`.

### OUTPUT (Required)

Variant calling output file. Can be VCF or BCF.

### REFERENCE (Required)

Reference genome fasta file. Note that gridss caches the entire reference genome in memory so ensure
that the working memory available to gridss comfortably exceeds the reference genome size.
16GB is recommended for multithreaded processing of whole-genome sequencing of mammalian genomes.

### INPUT (Required)

Input libraries. Specify multiple times (ie INPUT=file1.bam INPUT=file2.bam INPUT=file3.bam ) to process multiple libraries together.
Gridss considers all reads in each file to come from a single library.
Input files containing read groups from multiple different libraries should be split into an input file per-library.

### INPUT_READ_PAIR_MIN_CONCORDANT_FRAGMENT_SIZE, INPUT_READ_PAIR_MAX_CONCORDANT_FRAGMENT_SIZE

Per INPUT overrides of the default concordant fragment size calculation.
If your fragment size distribution is unusual, this option can be used to override the default value.

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

### WORKER_THREADS

Number of processing threads to use, including number of thread to use when invoking the aligner.
Note that the number of threads spawned by gridss is typically 2-3 times the number of worker threads due to asynchronous I/O threads
thus it is not uncommon to see over 100% CPU usage when WORKER_THREADS=1 as bam compression/decompression is a computationally expensive operation.
This parameter defaults to the number of cores available.

### WORKING_DIR

Directory to write intermediate results directories. By default, intermediate files for each input or output file are written to a subdirectory in the same directory as the relevant input or output file.
If WORKING_DIR is set, all intermediate results are written to subdirectories of the given directory.

### TMP_DIR

This field is a standard Picard tools argument and carries the usual meaning. Temporary files created during processes such as sort are written to this directory.

## libsswjni.so

Due to relatively poor performance of existing Java-based Smith-Waterman alignment packages, gridss incorporates a JNI wrapper to the striped Smith-Waterman alignment library [SSW](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library). lib/libsswjni.so is a precompiled linux version included in gridss If the precompiled version is not compatable with your linux distribution, the wrapper can be compiled by:

* git clone https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library
* make java

If libsswjni.so cannot be loaded, gridss will fall back to a slower java implementation. To ensure libsswjni.so is loaded, add the library using -Djava.library.path, or the LD_LIBRARY_PATH environment variable as per the JNI documentation.

### GRIDSS_CONFIG

Gridss uses a large number of configurable settings and thresholds which are, by default, hidden from the user.
These settings can be overriden by specifying a configuration file to use instead. Note that this configuration
file uses a different format to the Picard tools-compatable configuration file that is used instead of the
standard command-line arguments

When supplying a custom configuration, gridss will use the overriding settings for all properties specified
and fall back to the default for all properties that have not been overridden. Packaged within the gridss jar
is resources/gridss.properties (src/main/resource/gridss.properties in the gridss source code) which contains
a detailed listing of all configurable gridss settings.

# Output

Gridss is fundamentally a structural variation breakpoint caller. Variants are output as VCF breakends. Each call is a breakpoint consisting of two breakends, one from location A to location B, and a reciprocal record from location B back to A. Note that although each record fully defines the call, the VCF format required both breakend to be written as separate record.

## Quality score

Gridss calculates quality scores according to the model outlined in [paper].
As gridss does not yet perform multiple test correction or score recalibration, **QUAL scores are vastly overestimated for all variants**.
As a rule of thumb, variants with QUAL >= 1000 and have assembles from both sides of the breakpoint (AS > 0 & RAS > 0) are considered HIGH quality,
variant with QUAL >= 500 but can only be assembled from one breakend (AS > 0 | RAS > 0) are considered MEDIUM quality,
and variants with low QUAL score or lack any supporting assemblies are considered LOW quality. By default, the NO_ASSEMBLY filter is applied to all variants
without assembly support.

## Non-standard INFO fields

Gridss writes a number of non-standard VCF fields. These fields are described in the VCF header.

## BEDPE

Gridss supports conversion of VCF to BEDPE format using the VcfBreakendToBedpe utility program included in the gridss jar.

Calling VcfBreakendToBedpe with `INCLUDE_HEADER=true` will include a header containing column names in the bedpe file. These fields match the VCF INFO fields of the same name. For bedpe output, breakend information is not exported and per category totals (such as split read counts) are aggregrated to a single value.








