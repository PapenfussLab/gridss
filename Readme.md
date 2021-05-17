[![Release](https://img.shields.io/github/v/release/PapenfussLab/gridss)](https://github.com/PapenfussLab/gridss/releases)
[![Build Status](https://travis-ci.org/PapenfussLab/gridss.svg?branch=master)](https://travis-ci.org/PapenfussLab/gridss)
[![Language](http://img.shields.io/badge/language-java-brightgreen.svg)](https://www.java.com/)

# GRIDSS - the Genomic Rearrangement IDentification Software Suite

GRIDSS is a module software suite containing tools useful for the detection of genomic rearrangements. GRIDSS includes a genome-wide break-end assembler, as well as a structural variation caller for Illumina sequencing data. GRIDSS calls variants based on alignment-guided positional de Bruijn graph genome-wide break-end assembly, split read, and read pair evidence.

GRIDSS makes extensive use of the [standard tags defined by SAM specifications](http://samtools.github.io/hts-specs/SAMtags.pdf). Due to the modular design, any step (such as split read identification) can be replaced by another implementation that also outputs using the standard tags. It is hoped that GRIDSS can serve as an exemplar modular structural variant pipeline designed for interoperability with other tools.

If you have any trouble running GRIDSS, please raise an issue using the Issues tab above. Based on feedback from users, a user guide will be produced outlining common workflows, pitfalls, and use cases.

[Click here to download GRIDSS](https://github.com/PapenfussLab/gridss/releases)

Detailed documentation is being developed [here](https://github.com/PapenfussLab/gridss/wiki/GRIDSS-Documentation)

# Citation

For citing GRIDSS and for an overview of the GRIDSS algorithms, refer to our open access article: http://genome.cshlp.org/content/early/2017/11/02/gr.222109.117.abstract

Daniel L. Cameron, Jan Schröder, Jocelyn Sietsma Penington, Hongdo Do, Ramyar Molania, Alexander Dobrovic, Terence P. Speed and Anthony T. Papenfuss.
GRIDSS: sensitive and specific genomic rearrangement detection using positional de Bruijn graph assembly.
Genome Research, 2017
doi: 10.1101/gr.222109.117


For single breakend, structural variant phasing, and somatic variant calling, refer to our preprint:

GRIDSS2: harnessing the power of phasing and single breakends in somatic structural variant detection
https://www.biorxiv.org/content/10.1101/2020.07.09.196527v1

VIRUSBreakend usage should cite "VIRUSBreakend: Viral Integration Recognition Using Single Breakends" https://doi.org/10.1093/bioinformatics/btab343

# Pre-requisites

To run GRIDSS the following must be installed:

* java 1.8 or later
* R 4.0 or later
  * `gridss_somatic_filter` and `gridss_extract_overlapping_fragments` require the following R libraries:
    * argparser
    * tidyverse
    * stringdist
    * testthat
    * stringr
    * StructuralVariantAnnotation
    * rtracklayer
    * BSgenome package for your reference genome (optional)
* samtools 1.10 or later
* bwa
* bash
* getopt(1) (part of [util-linux](https://en.wikipedia.org/wiki/Util-linux))

To run VIRUSBreakend, kraken2, or repeatmasker annotations, the following additional software must be installed:
* kraken2
  * Note that `virusbreakend-build` requires all `kraken2-build` dependencies
* RepeatMasker
* bcftools

## Building GRIDSS

GRIDSS is mostly written in Java thus local building is not required.
Just [download the latest release](https://github.com/PapenfussLab/gridss/releases) and ensure you have the Pre-requistes installed.

If you wish to contribute to GRIDSS development, it can be built from source using maven with `mvn package`.

A prebuilt docker image is available as `gridss/gridss:latest` so building a docker image yourself is not necessary.
If you do wish to build the docker image yourself, you need to run `scripts/dev/create_release.sh` so the necessary build artifacts exist in your environment.

### Building gridsstools

Some performance-critical steps are implemented in C using htslib.
A precompiled version of `gridsstools` for linux is included as part of GRIDSS releases.
If this precompiled version does not run on your system you will need to build it from source.

To build `gridsstools` from source run the following:
```
git clone http://github.com/PapenfussLab/gridss/
cd gridss
git submodule init
git submodule update
cd src/main/c/gridsstools/htslib/
autoheader
autoconf
./configure && make
cd ..
autoheader
autoconf
./configure && make all
```

## Conda issues

Compiling with a conda environment active is likely to cause problems such as `undefined reference to 'libdeflate_crc32'`. This happens when the conda environment includes copies of the libraries used by htslib (z m bz2 lzma curl crypto pthread) without also including the headers for the libraries. This causes gridsstools to compile against the system headers, but link against the conda libraries, hence the error.

Run either `conda install htslib` or `conda deactivate` if you have problems compiling gridsstools in a conda environment.

# Running

Pre-compiled binaries are available at https://github.com/PapenfussLab/GRIDSS/releases. GRIDSS invokes external tools at multiple points during processing. By default this is bwa mem, but can be configured to use bowtie2 or another aligner.

The following programs are included in GRIDSS releases:

|program|description
|---|---|
gridss|GRIDSS assembler and structural variant caller. Use this to generate a GRIDSS SV VCF.
gridss_extract_overlapping_fragments|Extracts all alignments for read pairs with at least one aligment overlapping set of regions of interest. Correctly handles supplementary alignments. Use this script to extract reads of interest for targeted GRIDSS variant calling.
gridss_annotate_vcf_repeatmasker|Annotates breakpoint and single breakend inserted sequences with the RepeatMasker classification of the sequence.
gridss_annotate_vcf_kraken2|Annotates breakpoint and single breakend inserted sequences with the Kraken2 classification of the sequence.
virusbreakend|[See VIRUSBreakend README](https://github.com/PapenfussLab/gridss/blob/master/VIRUSBreakend_Readme.md)
virusbreakend-build|[See VIRUSBreakend README](https://github.com/PapenfussLab/gridss/blob/master/VIRUSBreakend_Readme.md)
gridss_somatic_filter|Somatic filtering script. Note that this has an equivalent java implementation in [GRIPSS](https://github.com/hartwigmedical/hmftools/tree/master/gripss).
gridsstools|C/htslib implementation of performance-critical steps. Currently used by `virusbreakend` and `gridss_extract_overlapping_fragments`

## gridss command-line arguments

```
Usage: gridss --reference <reference.fa> --output <output.vcf.gz> --assembly <assembly.bam> [--threads n] [--jar gridss.jar] [--workingdir <directory>] [--jvmheap 30g] [--blacklist <blacklist.bed>] [--steps All|PreProcess|Assemble|Call] [--configuration gridss.properties] [--maxcoverage 50000] [--labels input1,input2,...] input1.bam [input2.bam [...]]
```

argument|description
---|---
-o, --output|output VCF
-a, --assembly|location of the GRIDSS assembly BAM. This file will be created by GRIDSS
-r, --reference|reference genome to use. Must have a .fai index file and a bwa index
-t, --threads|number of threads to use. Defaults to the number of cores available.
-j, --jar|location of GRIDSS jar
-b/--blacklist|BED file containing regions to ignore. The ENCODE DAC blacklist is recommended for hg19. (Optional)
--jvmheap|size of JVM heap for assembly and variant calling. Defaults to 30g to ensure GRIDSS runs on cloud instances with 32gb memory.
--maxcoverage|maximum coverage. Regions with coverage in excess of this are ignored. (Default: 50000)
--labels|comma separated labels to use in the output VCF for the input files. Must have same number of entries as there are input files. Input files with the same label are aggregated (useful for multiple sequencing runs of the same sample). Labels default to input filenames, unless a single read group with a non-empty sample name exists in which case the read group sample name is used (which can be disabled by \"useReadGroupSampleNameCategoryLabel=false\" in the configuration file). If labels are specified, they must be specified for all input files.
--steps|processing steps to run. Defaults to all steps. Multiple steps are specified using comma separators. Available steps are preprocess,assemble,call. Useful to improve parallelisation on a cluster as preprocess of each input file is independent, and can be performed in parallel, and has lower memory requirements than the assembly step.
--jobindex|zero-based index of this assembly job node. Used to spread GRIDSS assembly across multiple compute nodes. Use only with `-s assemble`. Once all jobs have completed, a `-s assemble` or `-s all` job should be run to gather the results together.
--jobnodes|total number of assembly jobs scheduled.

The following additional optional arguments may be useful if GRIDSS fails to run in your environment, or you want to run with non-standard parameters.

argument|description
---|---
-c, --configuration|configuration file use to override default GRIDSS settings
--externalaligner|use the system version of bwa instead of the in-process version packaged with GRIDSS
--picardoptions|additional standard Picard command line options. Useful options include VALIDATION_STRINGENCY=LENIENT and COMPRESSION_LEVEL=0. See https://broadinstitute.github.io/picard/command-line-overview.html
--useproperpair|use SAM 'proper pair' flag to determine whether a read pair is discordant. Default: use library fragment size distribution to determine read pair concordance
--concordantreadpairdistribution|portion of read pairs distribution considered concordantly mapped. Default: 0.995
--keepTempFiles|keep intermediate files. Not recommended except for debugging due to the high disk usage.
--nojni|do not use JNI native code acceleration libraries (snappy, GKL, ssw, bwa).
	
_Warning_: all somatic R scripts treat the first bam file to be the matched normal, and any subsequent as tumour sample. If you are doing somatic calling, make sure you follow this convention.

### gridss steps

The following GRIDSS steps can be specified:
step|description
---|---
setupreference|Once-off setup generating additional files in the same directory as the reference. WARNING: multiple instances of GRIDSS attempting to perform `setupreference` at the same time will result in file corruption. Make sure these files are generated before running parallel GRIDSS jobs.
preprocess|Pre-process input BAM files. Can be run per input file.
assemble|Perform GRIDSS breakend assembly. Can split up across multiple nodes using `--jobindex` and `--jobnodes`.
call|Perform variant calling.
all| Run all steps (Default)

At present, command line valiation is performed independently of which steps are run. When splitting GRIDSS into multiple cluster jobs, the same command line parameters should be specified for every job except for:
- `--input` on preprocess jobs (one input per file)
- `--jobindex` and `--jobnodes` on assembly jobs

# FAQ

### How do I run GRIDSS on multiple samples?

Just specify multiple BAMs on the command line. GRIDSS will perform joint calling and provide a per-BAM breakdown of support.

### Should I do joint calling or run each sample individually?

**Joint calling should always be used for related samples** (e.g. tumour/normal or trio calling).
Joint calling will ensure that a common variant near the single-sample threshold of detection will be reliably reported as a shared variant.
This is not the case if the calling were done individually.
Note that this particular behaviour is not specific to GRIDSS and is common to all variant callers (hence the joint calling support in many callers).

Joint calling allows for sensivity detection of variants that are present subclonally (or at low coverage) that would not be detected if called individually.

GRIDSS performs joint assembly then reports a per-sample breakdown. Joint calling has higher coverage of shared variants thus resulting in more reliable assembly of that variant.

Determining whether two SV calls in two different VCFs are actually the same call is non-trivial.
Imprecise calls are especially problematic since the coordinates may differ between the VCF, or a call may be precise in one VCF and not in the other.
A good example of why reconciling SV calls is so problematic is the case where call A (chrX:1-99->chrY:1-99) overlaps call B (chrX:50-149->chrY:1-99), call B overlaps call C (chrX:100-199->chrY:1-99), but A does not overlap C at all. Joint calling obviates this step.

### How do I perform tumour/normal somatic variant calling?

Jointly call on all samples from the patient.
It is strongly recommended that the normal be the first argument as that is what downstream steps expect.
For example,  `gridss ... patientX_normal.bam patientX_primary.bam patientX_met.bam`.
To filter to somatic calls, use the `gridss_somatic_filter` script included in the GRIDSS release.

### What aligner should I use?

The default of `bwa mem` is sufficient for most use cases.

Although GRIDSS aims to be aligner agnostic, not all aligners output BAM files suitable for processing by GRIDSS. GRIDSS requires:
* One alignment per read. Supplementary (split read) alignments are ok, but secondary alignments are not.
  * This means that aligner settings such as the `-a` option of bwa mem and the `-k` and `-a` options of bowtie2 are unsuitable.
* MAPQ to meaningfully follow the SAM specifications. Aligners that do not follow the specifications (e.g. subread) will have worse results.

Options such as the `-Y` option of bwa mem, or the fact that bowtie2 does not do split read alignment are not problematic as these differences are corrected in the GRIDSS preprocessing step.

### How do I tell GRIDSS multiple BAMs are from the same sample?

Use the `--labels` command line option. Eg: `--labels sample1,sample1,sample2 sample1_library1.bam sample1_library2.bam sample2.bam`

### Why are there ALT alleles with `.` in the output?

This is the [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) notation for single breakend variant calls. See section 5.4.9 of the specifications document. These calls indicate that a breakpoint was found at this location but the partner location could not be unambiguously determined.

### How do I get output like the GRIDSS PURPLE LINX figures?

Run the [integrated GRIDSS PURPLE LINX pipeline script](https://github.com/hartwigmedical/gridss-purple-linx) or the docker image gridss/gridss-purple-linx:latest.

### How do I do RepeatMasker annotation of breakend sequences?

Run `gridss_annotate_vcf_repeatmasker` on the GRIDSS output.

### How do I do viral annotation?

Use VIRUSBreakend for viral annotations. See the [VIRUSBreakend README](https://github.com/PapenfussLab/gridss/blob/master/VIRUSBreakend_Readme.md) for more details.

### What does `gridss_somatic_filter` actually do?

See documentation at https://github.com/PapenfussLab/gridss/wiki/Somatic-Filtering

The Hartwig Medical Foundation has reimplemented `gridss_somatic_filter` in Java as [GRIPSS](https://github.com/hartwigmedical/hmftools/tree/master/gripss).
GRIPSS is much faster, has additional features, and is the recommended tool for somatic filtering of GRIDSS output.

### How do I create the panel of normals required by `gridss_somatic_filter`?

If you are using hg19 or hg38, then a PON based on Dutch samples is available from https://resources.hartwigmedicalfoundation.nl/.
Make sure the reference you are using and the PON both use the same chromosome notation or nothing will get filtered (e.g. `1` vs `chr1`).
If these are not appropriate, you'll need to create your own using `gridss.GeneratePonBedpe`

Here is an example that generates a PON from every VCF in the current directory:

```
mkdir -p pondir
java -Xmx8g \
	-cp ~/dev/gridss/target/gridss-2.10.2-gridss-jar-with-dependencies.jar \
	gridss.GeneratePonBedpe \
	$(ls -1 *.vcf.gz | awk ' { print "INPUT=" $0 }' | head -$n) \
	O=pondir/gridss_pon_breakpoint.bedpe \
	SBO=pondir/gridss_pon_single_breakend.bed \
	REFERENCE_SEQUENCE=$ref_genome
```

The score fields of the bedpe/bed files is the count of the number of samples that variant was found in.
I recommended filtering these output files to only variants found in 3+ samples.

Note that `gridss_somatic_filter` requires the files to be named `gridss_pon_breakpoint.bedpe` and `gridss_pon_single_breakend.bed`.

### How do I merge PONs?

Merging of PONs is not supported but incremental updates are.
Use `INPUT_BEDPE` and `INPUT_BED` files to `gridss.GeneratePonBedpe`.
For example, you can add your samples to the Hartwig PONs by pointing `INPUT_BEDPE` and `INPUT_BED` to the Hartwig PONs and adding your VCFs as outlined above.

### Should I process each input BAM separately or together?

Wherever possible, samples should be processed together. Joint calling enables the detection of low allelic fraction SVs, as well as making the downstream analysis much easier (a single VCF with a breakdown of support per sample is much easier to deal with than multiple VCFs - the matching logic required to determine if two SVs are equivalent is non-trivial).

GRIDSS joint calling has been tested on up 12 samples with ~1000x aggregate coverage. If you have hundreds of samples, joint assembly will likely be computationally prohibitive and you will need to perform assembly in batches, then merge the results together.

WARNING: multiple instances of GRIDSS generating reference files at the same time will result in file corruption. Make sure `setupreference` files have been generated before running parallel GRIDSS jobs.

#### How do I perform assembly on multiple nodes?

GRIDSS preprocessing and assembly can be spread across multiple nodes.
Preprocess parallelisation is one job per input file.
Assembly parallelisation is one thread per genomic region.
To reduce wall times, regions can be distributed across multiple nodes using the `--jobindex` and `--jobnodes` parameters.

Here is an example:
```
gridss -s setupreference # once-off-setup
# in parallel:
gridss -s preprocess input1.bam
gridss -s preprocess input2.bam
gridss -s preprocess input3.bam
gridss -s preprocess input4.bam
# wait for all preprocessing jobs to complete
# Perform assembly parallel (across three nodes in this example)
gridss -s assemble --jobindex 0 --jobnodes 3 -a assembly.bam  input1.bam input2.bam input3.bam input4.bam
gridss -s assemble --jobindex 1 --jobnodes 3 -a assembly.bam input1.bam input2.bam input3.bam input4.bam
gridss -s assemble --jobindex 2 --jobnodes 3 -a assembly.bam input1.bam input2.bam input3.bam input4.bam
# wait for all assembly jobs to complete
# Gather the assembly results togther. This job is essentially a file copy and completes very quicky
gridss -s assemble -a assembly.bam input1.bam input2.bam input3.bam input4.bam
# perform variant calling
gridss -s call  -a assembly.bam input1.bam input1.bam input2.bam input3.bam input4.bam
```

####  How do I perform batched assembly?

If you can perform joint assembly, do so.
If you run out of memory when performing joint assembly, or hit regularly encounter assembly timeouts, batched assembly is required.
You are likely to encounter these issues at, or above, 1000x coverage.

To perform batched assembly, run the GRIDSS `assemble` step multiple times each with only a subset of the input files
If you use input labels, all inputs with the same label must be in the same batch.

Here is an example:
```
gridss -s setupreference # once-off-setup
gridss -s preprocess input1.bam
gridss -s preprocess input2.bam
gridss -s preprocess input3.bam
gridss -s preprocess input4.bam
gridss -s assemble -a assembly12.bam input1.bam input2.bam
gridss -s assemble -a assembly34.bam input3.bam input4.bam
gridss -s call -a assembly12.bam -a assembly34.bam input1.bam input2.bam input3.bam input4.bam
```

Related samples should always be assembled together in the same batch.

### I encountered an error. What should I do?

* Check the bottom of this page for commonly encountered errors and their solutions

### How many threads should I use?

GRIDSS has been optimised to run on a 8core/32gb cloud compute node.

If scaling above 8 cores, it is recommended to run multiple GRIDSS assembly processes and use the `--jobindex` and `--jobnodes` parameters. with each job allocated 8 cores/32gb.

### How much memory should I give GRIDSS?

GRIDSS has been optimised to run on a 8core/32gb cloud compute node.

At least 4GB + 2GB per thread. It is recommended to run GRIDSS with max heap memory (-Xmx) of 8GB for single-threaded operation
(WORKER_THREADS=1), 16GB for multi-core desktop operation, and 31GB for server operation. Note that due to Java's use of [Compressed Oops](http://docs.oracle.com/javase/7/docs/technotes/guides/vm/performance-enhancements-7.html#compressedOop), specifying a max heap size of between 32-48GB effectively reduces the memory available to GRIDSS so is strongly discouraged.

### Why does GRIDSS use more CPU than the limit specified with `--threads`?

GRIDSS has been optimised to run on a 8core/32gb cloud compute node.

`--threads` specifies the size of the worker thread pool. IO, BAM decompression, and parsing are in their own thread pool which is not part of the worker thread pool. The pre-processing, variant calling, and annotation steps also perform some work that is executed in dedicated threads independent of the worker thread pool. Combined, this approach means that max CPU utilisation can exceed the thread count specified.

Asynchronous IO defaults can be changed by editing the `jvm_args` argument in `gridss`.

### Should I include alt contigs in the reference?

GRIDSS relies on the aligner to determine the mapping location and quality of assembly contigs and to identify split read from soft clipped reads. GRIDSS considers alignments with low mapping quality (default mapq <= 10) to not be uniquely aligned and treats them as unaligned. If alt contigs are included in an aligner that is not alt-aware then hits to sequences that are in both the primary reference contigs and the alt contigs will be given a low mapq by the aligner and will be treated as unaligned by GRIDSS. By default, GRIDSS uses bwa mem for alignment so by including the bwa alt contig definition file, reads from regions with alt homology will be preferentially aligned to the reference with a correspondingly improved mapq. Whether or not to include alt contigs depends on what sort of downstream analysis you intend to perform and how your intend to handle structural variants involving alt contigs. That said, if your reference genome includes alt contigs and a bwa alt contig definition file is available for your genome, you should use it.

### How do I process only my region of interest?

Extract all fragments overlapping your region of interest, using `gridss_extract_overlapping_fragments` then run `gridss` on the subset bam.
`gridss_extract_overlapping_fragments` is almost identical to filtering using `samtools view` except that it extracts alignment records for any fragment overlapping a region of interest (ie mate reads and supplementary alignments).
That is, all records with read names matching the read name of an alignment overlapping any region of interest.

### How do I use GRIDSS to validate the calls from another caller?

- Convert your calls to a BED file containing the start and end positions of all SVs called by your other caller.

- Expand intervals by at least 10kbp. Too small a window will have a negative impact on GRIDSS QUAL scores (since they're emperically weighted, taking only regions with soft clipped reads will cause GRIDSS to massively downweight soft clips when scoring). Alternatively, for unbiased scoring, run `gridss.CollectGridssMetrics` on your input file and rename the `.gridss.working` directory to the name of your targeted bam file to enable the targeted bam to use the full bam metrics.

- Process as per the region of interest processing outlined above.

### Why does GRIDSS use centre-aligment?

'Normalised' SNV/indel calls are left-aligned. Why does GRIDSS not use this convention? There are two reasons:

For ++ or -- breakpoints, left-aligning the lower breakend will force right-alignment of the upper breakend. Similarly for right-alignment of the lower. This means that it is impossible to univerally left-align, or right-align breakpoints without resulting in an incorrect nominal call position. Centre-alignment is the option that does not cause this problem ( technically speaking, you still have an off-by-one problem for odd-length homology but that's less problematic).

Secondly, using left or right alignment for imprecise call will result in the nominal call being at the edge of the confidence interval bounds. Centre-aligning imprecise calls makes sense as it is (usually) the centre position that is the most likely to be correct.

### What does `gridss_annotate_kraken2` output?

`gridss_annotate_kraken2` adds Kraken2 classifications to single breakend and breakpoint inserted sequences.
The [NCBI taxonomy ID](https://www.ncbi.nlm.nih.gov/taxonomy) for the inserted sequences is in the `INSTAXID` INFO field.

## GRIDSS JAR

GRIDSS takes a modular approach and the GRIDSS jar consists of a collection of separate tools. Each tool in the GRIDSS pipeline can be run independently. The following data flow diagram gives an overview of the GRIDSS pipeline used when running `gridss`.

![GRIDSS data flow diagram](https://docs.google.com/drawings/d/1aXFBH0E9zmW4qztHIEliZfsLCHJa6_-l624Frq1X-Ms/pub?w=973&h=760)

#### CallVariants

This tool runs every step in the variant calling pipeline. This entry point has been superceeded by `gridss` but is retained for backward compatibility with existing pipeline. `gridss` is preferred as it has lower peak memory usage and is slightly faster due to the use of samtools for sorting instead of htsjdk.

#### CollectGridssMetrics

GRIDSS requires a number of input library metrics to be calculated. The metrics calculations programs are invoked in the same manner as Picard tools metrics. This program functions similarly to Picard tools CollectMultipleMetrics but, by default, only extracts the metrics required by GRIDSS.

#### CollectCigarMetrics

Collects metrics regarding the size and type of read alignment CIGAR elements.

#### CollectIdsvMetrics

Collects generic library-level metrics such as the read and read pair mapping rates.

#### CollectMapqMetrics

Collects metrics regarding distribution of read mapping quality scores.

#### CollectTagMetrics

Collects metrics regarding presence of SAM tags.

#### ExtractSVReads

Extracts the subset of reads that provide potential support for structural variation. These reads fall into one or more of the following classes:
* Reads aligned with an insertion or deletion in the read alignment
* Soft clipped reads
* Split reads (identified by the SA SAM Tag)
* Discordant read pairs
* Read pairs with only 1 read mapped
* Unmapped reads

#### SoftClipsToSplitReads

Identifies split reads by iterative realignment of soft clipped bases with a NGS aligner. By default, bwa mem is used for read alignment.
The GRIDSS pipeline runs this on all input files, as well as the breakend assembly contigs generated by AssembleBreakends.

#### ComputeSamTags

This step recomputes redundant information to ensure data consistency and writes additional standards-defined SAM tags by ppulates SAM tags required in downstream steps, softening hard clips, and correcting data inconsistencies in the input files. In particular, `AssembleBreakends` requires the `R2` tag to be populated. Since `R2` is populated from the mate sequence and softening hard clips requires all split read alignment records, this step requires all records with the same read name to be adjacent in the BAM file (ie read name grouped sort order). Data inconsistencies corrected include:

- `NM` tag not matching edit distance to reference
- Mate CIGAR and alignment field not matching actual alignment and CIGAR of mate (e.g. GATK indel realignment does not update mate information when realigning reads)
- `SA` tag not matching supplementary alignment
- Supplementary alignments not duplicated marked (markdup only marks primary alignments)
- Supplementary alignments incorrectly flagged as secondary alignments

This step can recalculate the following SAM tags:

- SA (default)
- NM (default)
- R2 (default)
- Q2
- MC (default)
- MQ (default)
- CC
- CP
- HI
- IH
- FI

#### AssembleBreakends

Generates breakend assemblies from the input reads. Breakend assemblies are written as synthetic soft clipped reads.
Each of these reads corresponds to a breakend assembly contig composed of reads anchored (either directly, or via the mate read)
to the location of the breakend assembly.

Breakend contigs composed entirely of discordant read pairs or reads with unmapped mates cannot be uniquely placed as there exists
an interval over which the contig is anchored. These alignments are written using a placeholder CIGAR alignment of the form XNX.
For example, a breakend contig read with an alignment of 1X50N1X150S represents a 150bp contig in which the breakend is expected to
occur at one of the 52 genomic positions given by the placeholder XNX alignment interval.

Breakpoints are identified by running SoftClipsToSplitReads on the breakend assembly contigs. High quality breakpoints are expected
to have two independent breakend assemblies (one from each breakend of the breakpoint).

#### IdentifyVariants

Identifies putative structural variants from the reads providing potential SV support, and the breakend assembly contigs.

#### AnnotateVariants

Annotates breakpoint calls performing AllocateEvidence, AnnotateReferenceCoverage, AnnotateInexactHomology annotation.

##### AllocateEvidence

Uniquely allocate reads/read pairs to identified breakpoints by ensuring that for each read/read pair, only the single
best alignment is retained (relevant for input files containing multiple mapping locations for each read).
Read counts and other summary statistics are collated.

##### AnnotateReferenceCoverage

Calculates the number of reads and read pairs providing support for the reference allele at each breakpoint.

#### AnnotateInexactHomology

Calculates the size of the inexact homology between the reference sequence and the breakpoint sequence. Breakpoints
with long inexact homology are possibly due to alignment artifacts causing false positive breakpoint calls between
regions of homologous sequence.

#### AnnotateInsertedSequence

Finds potential mapping locations for single breakends and breakpoint insert sequences in the given reference. Used for RepeatMasker annotation, and viral integration detection.

#### GeneratePonBedpe

This tool aggregates variants across multiple VCFs and counts the number of samples supporting each variant. Only the first sample per VCF is processed which is useful for generating a panel of normals (PON) from a cohort of cancer samples with matched normals. Output is a bedpe (breakpoint) and bed (single breakend) file suitable for use by `gridss_somatic_filter`.

## Common Parameters

GRIDSS programs have a large number of parameters that can be be adjusted. The default parameter set has been tested with paired-end Illumina data ranging from 2x36bp to 2x250bp and should give a reasonably good result. Command line used parameters are listed below.

### OUTPUT (Required)

Variant calling output file. Can be VCF or BCF.

### REFERENCE_SEQUENCE (Required)

Reference genome FASTA file. GRIDSS requires that the reference genome supplied exactly matches
the reference genome of all input files.
The reference genome must be in FASTA format and must have a tabix (.fai) index and an
index for the NGS aligner (by default bwa). The NGS aligner index prefix must match
the reference genome filename. For example, using the default setting against the reference
file reference.fa, the following files must be present and readable:

File | Description
------- | ---------
reference.fa | reference genome
reference.fa.fai | Tabix index
reference.fa.amb | bwa index
reference.fa.ann | bwa index
reference.fa.bwt | bwa index
reference.fa.pac | bwa index
reference.fa.sa | bwa index

These can be created using `samtools faidx reference.fa` and `bwa index reference.fa`

A .dict sequence dictionary is also required but GRIDSS will automatically create one if not found.

### INPUT (Required)

Input libraries. Specify multiple times (i.e. `INPUT=file1.bam INPUT=file2.bam INPUT=file3.bam`) to process multiple libraries together.

Input files must be coordinate sorted SAM/BAM/CRAM files.

GRIDSS considers all reads in each file to come from a single library.

Input files containing read groups from multiple different libraries should be split into an input file per-library.

The reference genome used for all input files should match the reference genome supplied to GRIDSS.

### INPUT_LABEL

Labels to allocate inputs. The default label for each input file corresponds to the file name but can be overridden by
specifying an INPUT_LABEL for each INPUT. The output for any INPUT files with the same INPUT_LABEL will be merged.

### ASSEMBLY (Required)

File to write breakend assemblies to. It is strongly recommended that the assembly filename corresponds to the OUTPUT filename. Using ASSEMBLY=assembly.bam is problematic as (like the INPUT files) the assembly file is relative not to WORKING_DIR, but to the current directory of the calling process. This is likely to result in data corruption when the same assembly file name is used on different data sets (for example, writing assembly.bam to your home directory when running on a cluster).

### BLACKLIST

BED blacklist of regions to exclude from analysis. The [ENCODE DAC blacklist](https://www.encodeproject.org/annotations/ENCSR636HFF/)
is recommended when aligning against hg19.

Unlike haplotype assemblers such as TIGRA and GATK, GRIDSS does not abort assembly when complex assembly graphs are encountered. Processing of these graphs slows down the assembly process considerably, so if regions such as telomeric and centromeric regions are to be excluded from downstream analysis anyway, assembly of these regions is not required. It is recommended that a blacklist such as the [ENCODE DAC blacklist](https://www.encodeproject.org/annotations/ENCSR636HFF/) be used to filter such regions. Inclusion of additional mappability-based blacklists is not required as GRIDSS already considers the read mapping quality.

### READ_PAIR_CONCORDANT_PERCENT

Portion (0.0-1.0) of read pairs to be considered concordant. Concordant read pairs are considered to provide no support for structural variation.
Clearing this value will cause GRIDSS to use the 0x02 proper pair SAM flag written by the aligner to determine concordant pairing.
Note that some aligners set this flag in a manner inappropriate for SV calling and set the flag for all reads with the expected orientation and strand regardless of the inferred fragment size.

### INPUT_MIN_FRAGMENT_SIZE, INPUT_MAX_FRAGMENT_SIZE

Per input overrides for explicitly specifying fragment size interval to be considered concordant. As with INPUT_LABEL, these must be specified
for all input files. Use null to indicate an override is not required for a particular input (e.g.
`INPUT=autocalc.bam INPUT_MIN_FRAGMENT_SIZE=null INPUT_MAX_FRAGMENT_SIZE=null INPUT=manual.bam INPUT_MIN_FRAGMENT_SIZE=100 INPUT_MAX_FRAGMENT_SIZE=300` )

### WORKER_THREADS

Number of processing threads to use, including number of threads to use when invoking the aligner.
Note that the number of threads spawned by GRIDSS is greater than the number of worker threads due to asynchronous I/O threads thus it is not uncommon to see over 100% CPU usage when WORKER_THREADS=1 as bam compression/decompression is a computationally expensive operation.
This parameter defaults to the number of cores available.

### WORKING_DIR

Directory to write intermediate results directories. By default, intermediate files for each input or output file are written to a subdirectory in the same directory as the relevant input or output file.
If WORKING_DIR is set, all intermediate results are written to subdirectories of the given directory.

### TMP_DIR

This field is a standard Picard tools argument and carries the usual meaning. Temporary files created during processes such as sort are written to this directory.

### samjdk defines

GRIDSS uses [htsjdk](https://github.com/samtools/htsjdk) as a SAM/BAM/CRAM/VCF parsing library. The following htsjdk java command-line options are strongly recommended for improved performance:

* -Dsamjdk.use_async_io_read_samtools=true
* -Dsamjdk.use_async_io_write_samtools=true
* -Dsamjdk.use_async_io_write_tribble=true

## libsswjni.so

Due to relatively poor performance of existing Java-based Smith-Waterman alignment packages, GRIDSS incorporates a JNI wrapper to the striped Smith-Waterman alignment library [SSW](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library). GRIDSS will attempt to load a precompiled version which is supplied as part of the GRIDSS package (a libsswjni.so file will be created in the TMP_DIR when GRIDSS is run). If the precompiled version is not compatible with your linux distribution, or you are running a different operating system, recompilation of the wrapper from source will be required. When recompiling, ensure the correct libsswjni.so is loaded using -Djava.library.path, or the LD_LIBRARY_PATH environment variable as per the JNI documentation.

If you have an older CPU that does not support SSE instructions, GRIDSS will terminate with a fatal error when loading the library. Library loading can be disabled by adding `-Dsswjni.disable=true` to the GRIDSS command line. If libsswjni.so cannot be loaded, GRIDSS will fall back to a (50x) slower java implementation which will result in the GRIDSS inexact homology variant annotation step running very slowly.

### CONFIGURATION_FILE

GRIDSS uses a large number of configurable settings and thresholds which for ease of use are not included
as command line arguments. Any of these individual settings can be overriden by specifying a configuration
file to use instead. Note that this configuration file uses a different format to the Picard tools-compatible
configuration file that is used instead of the standard command-line arguments.

When supplying a custom configuration, GRIDSS will use the overriding settings for all properties specified
and fall back to the default for all properties that have not been overridden. Details on the meaning
of each parameter can be found in the javadoc documentation of the `au.edu.wehi.idsv.configuration` classes.

# Output

GRIDSS is fundamentally a structural variation breakpoint caller. Variants are output as VCF breakends. Each call is a breakpoint consisting of two breakends, one from location A to location B, and a reciprocal record from location B back to A. Note that although each record fully defines the call, the VCF format requires both breakends to be written as separate records.

To assist in downstream analysis, the `StructuralVariantAnnotation` R BioConductor package is strongly recommended.
Operations such as variant filter, annotation and exporting to other formats such as BEDPE can be easily accomplished using this package in conjuction with the BioConductor annotation packages.

## Quality score

GRIDSS calculates quality scores according to the model outlined in the [paper](http://biorxiv.org/content/early/2017/02/21/110387).
As GRIDSS does not yet perform multiple test correction or score recalibration, QUAL scores are vastly overestimated for all variants.
As a rule of thumb, variants that have QUAL >= 1000 and have assemblies from both sides of the breakpoint (AS > 0 & RAS > 0) are considered of high quality,
variants with QUAL >= 500 but that can only be assembled from one breakend (AS > 0 | RAS > 0) are considered of intermediate quality,
and variants with low QUAL score or lack any supporting assemblies are considered to be of low quality.

## Non-standard INFO fields

GRIDSS writes a number of non-standard VCF fields. These fields are described in the VCF header.

## BED/BEDPE export

The recommended way to convert GRIDSS output to BEDPE is via the BioConductor `StructuralVariantAnntotation` package.

```
library(StructuralVariantAnnotation)
library(rtracklayer)

vcf = readVcf("gridss.vcf")

# Export breakpoints to BEDPE
bpgr = breakpointRanges(vcf)
# TODO: add your event filtering here. The default GRIDSS output is very verbose/sensitive.
write.table(breakpointgr2bedpe(bpgr), file="gridss_breakpoints.bedpe", sep="\t", quote=FALSE, col.names=FALSE)
	
# Export single breakends to BED
begr = breakendRanges(vcf)
# TODO: add your event filtering here. The default GRIDSS output is very verbose/sensitive.
begr$score = begr$QUAL
export(begr, con="gridss_single_breakends.bed")

```

## Visualisation of results

When performing downstream analysis on variant calls, it can be immensely useful to be able to inspect
the reads that the variant caller used to make the variant calls. As part of the GRIDSS pipeline, the following
intermediate files are generated:

* *INPUT*.sv.bam
* *ASSEMBLY*.sv.bam

The inputsv  file contains all reads GRIDSS considered as providing putative support for
any potential breakpoint, including breakpoints of such low quality that GRIDSS did not
make any call. This file includes all soft clipped, indel-containing, and split reads, as
well as all discordant read pairs and pairs with only one read mapped.

Split reads can be identified by the presence of a
[SA SAM tag](http://samtools.github.io/hts-specs/SAMtags.pdf).

GRIDSS treats breakend assemblies as synthetic soft clipped read alignments thus assemblies
are displayed in the same manner as soft clipped/split reads.

* No SA tags indicates a breakpoint could not be unambiguously identified from the breakend contig.
* The source of the breakend contig can be identified by the read without the "Supplementary alignment" 0x800 SAM flag set.
* More than two alignments for a breakend contig indicates the contig spans a complex event involving multiple breakpoints.

## Intermediate Files

GRIDSS writes a large number of intermediate files. If rerunning GRIDSS with different parameters on the same input, all intermediate files must be deleted, or a different WORKING_DIR specified. All intermediate files are written to the WORKING_DIR directory tree, with the exception of temporary sort buffers which are written to TMP_DIR and automatically deleted at the conclusion of the sort operation.

File | Description
------- | ---------
gridss.* | Temporary intermediate file
gridss.lock.breakend.vcf | Lock directory to ensure that only one instance of GRIDSS is running for any given output file.
*.bai | BAM index for coordinate sorted intermediate BAM file.
WORKING_DIR/*file*.gridss.working | Working directory for intermediate files related to the given file.
WORKING_DIR/*input*.gridss.working/*input*\*_metrics | Various summary metrics for the given input file.
WORKING_DIR/*input*.gridss.working/*input*.realign.*N*.fq | Split read identification fastq file requiring alignment by NGS aligner.
WORKING_DIR/*input*.gridss.working/*input*.realign.*N*.bam | Result of NGS alignment.
WORKING_DIR/*input*.gridss.working/*input*.sv.bam | Subset of input reads considered by GRIDSS. This file is useful for visualisation of the supporting reads contributing to structural variant calls. Note that this file includes split read alignments identified from soft clipped reads in the input file.
WORKING_DIR/*assembly*.gridss.working/*assembly*.sv.bam | Assembly contigs represented as soft clipped or split reads. These are the assembly contigs GRIDSS uses for variant calling. Note that, as with split reads, GRIDSS uses the SA tag to encode split read alignments.
WORKING_DIR/*output*.gridss.working/*output*.breakpoint.vcf | Raw unannotated variant calls before unique allocation of multi-mapping reads.

## Building from source

Maven is used for build and dependency management which simplifies compile to the following steps:

* `git clone https://github.com/PapenfussLab/gridss`
* `cd gridss`
* `mvn clean package`

If GRIDSS was built successfully, a combined jar containing GRIDSS and all required libraries located at target/GRIDSS-_VERSION_-gridss-jar-with-dependencies.jar will have been created.

# Error Messages

For some error messages, it is difficult to determine the root cause and what to do to fix it.
Here is a list of key phrases of errors encountered by users and their solution

### htsjdk.samtools.SAMFormatException: SAM validation error

Your input file does not conform to the SAM/BAM specifications. Solutions are to fix the input file so it conforms to the specifications (recommended) or add `--picardoptions VALIDATION_STRINGENCY=LENIENT` to ignore the error. Note that not all errors can be ignored.

### Aborting since lock gridss.lock._OUTPUT_ already exists. GRIDSS does not support multiple simulatenous instances running on the same data.

Multiple instances of GRIDSS were run on the same data. GRIDSS does not yet support MPI parallelisation across multiple machines. Use the WORKER_THREADS parameter to specify the desired level of multi-threading. If using a cluster/job queuing system, a single non-MPI job should be submitted and either WORKER_THREADS explicitly set to the number of cores associated with the job requests, or the job should request the entire node.

If the lock directory exists and you know a GRIDSS process is not running (eg: the GRIDSS process was killed), then you can safely delete the lock directory.

### Exception in thread "main" java.lang.UnsupportedClassVersionError: au/edu/wehi/idsv/Idsv : Unsupported major.minor version 52.0

You are attempting to run GRIDSS with an old Java version. GRIDSS requires Java 8 or later.

### ExternalProcessFastqAligner	Subprocess terminated with with exit status 1. Alignment failed for _INPUT_.realign.0.fq

The external aligner (bwa) could not be run. The most common causes of this are:
- bwa is not on `PATH`
 - Does running "bwa" print out the bwa usage message? If you are using a cluster, you may have to add bwa to your `PATH` (eg `module add bwa`).
- bwa index does not exist
- bwa index has incorrect suffix
 - e.g. if the reference is ref.fa the index must be ref.fa.bwt _not_ ref.bwt

Can you run the bwa command exactly as it appears in the error message?

###  (Too many open files)

GRIDSS has attempted to open too many files at once and the OS file handle limit has been reached.
On linux 'ulimit -n' displays your current limit. This error likely to be encountered if you have specified a large number of input files or threads. The following solution is recommended:
* Increase your OS limit on open file handles (eg `ulimit -n _<larger number>_`)
  * Note that many linux systems have a default hard limit on open file handles of 4096 which with many samples is frequently too still too few. Increasing the hard limit requires root access.
* Increase the chunk size. The default chunk size is 10 million bases. This can be increased by adding a `chunkSize=50000000` line a `gridss.properties` file and adding `CONFIGURATION_FILE=gridss.properties` to the GRIDSS command line. Note that this will increase the number of bases processed by each job thus reduce the level of parallelisation possible.
* Reduce number of worker threads. A large number of input files being processed in parallel results in a large number of files open at the same time.

If those options fail, your remaining options are:
* Added `-Dgridss.defensiveGC=true` to the java command-line used for GRIDSS. Memory mapped file handles are not released to the OS until the buffer is garbage collected . This option add a request for garbage collection whenever a file handle is no longer used. This is a significant overhead and is not a good option for sparse data samples (such as exome or targetted sequencing) - increasing the chunk size is a much better option for these samples.
* As a last-ditch effort, you can keep rerunning GRIDSS until it completes. If you are using the default entry point of `gridss.CallVariants` and have `-Dgridss.gridss.output_to_temp_file=true`, then you can rerun GRIDSS and it will continue from where it left off. Assuming it doesn't keep dying at the same spot, it will eventually complete.

### Reference genome used by _input.bam_ does not match reference genome _reference.fa_. The reference supplied must match the reference used for every input.

The reference genome used to align input.bam does not match the reference genome supplied to GRIDSS.
If the differences are purely based on chromosome name and ordering, the Picard tools utility ReorderBam
can be used to fix chromosome orderings.

### Unable to use sswjni library - assembly will be very slow. Please ensure libsswjni for your OS and architecture can be found on java.library.path

The sswjni library could not be loaded as the precompiled version is not compatable with your environment. See the sswjni sections for details on how to disable libsswjni or recompile it for your system.

### "Segmentation Fault", fatal JVM error, or no error message.

This is likely to be caused by a crash during alignment in libsswjni. See the sswjni sections for details on how to disable libsswjni or recompile it for your system.

### Illegal Instruction

Your CPU does not support the SSE2 instruction set. See the sswjni sections for details on how to disable libsswjni.

### Java HotSpot(TM) 64-Bit Server VM warning: INFO: os::commit_memory(0x00007fc36e200000, 48234496, 0) failed; error='Cannot allocate memory' (errno=12)

GRIDSS has run out of memory. Either not enough memory has been allocated to run GRIDSS or GRIDSS has attempted to memory map too many files (See "(Too many open files)"). In both cases, restart GRIDSS (increasing the memory available if required) and GRIDSS will continue from where it left off.

### java.lang.AssertionError: java.lang.ClassNotFoundException: com.sun.tools.javac.api.JavacTool

You are running GRIDSS in multi-mapping mode using only a JRE instead of a full JDK. Update your PATH and JAVA_HOME to a Java 1.8+ JDK installation.

### htsjdk.samtools.util.RuntimeIOException: java.io.IOException: No space left on device

Just like Picard tools and htsjdk libraries that GRIDSS uses, intermediate files are sorted according the the `TEMP` file location. On many system, /tmp does not have enough space to sort a BAM file so it is possible to run out of intermediate file storage even if you have plenty of space left on the file system the input and output files are stored on. Using the same command-line options as Picard tools, the intermediate files location can be set using the `TMP_DIR` command-line argument.

It's also possible that you've just run out of space.

# Contributing to GRIDSS

Bug fix pull requests are always welcome.

If you have a feature you would like to implement, please first raise an issue outlining the problem and proposed solution so as to avoid any wasted or duplicated effort.
To build GRIDSS, the development environment requires maven and all GRIDSS and VIRUSBreakend dependencies (see GRIDSS Dockerfile and VIRUSBreakend readme).
To run all test cases, a `../ref/` directory is also required containing the following genomes (with all associated GRIDSS `--setupreference` files created):

 - hg19.fa
 - Homo_sapiens_assembly38.fasta
 - hg38.fa
 - Escherichia_coli_bl21_de3_.ASM956v1.dna.toplevel.fa

`scripts/dev/` contains helpful utilities for GRIDSS development work.


# Acknowledgement

GRIDSS uses the YourKit Java Profiler for Performance Tuning.

![alt text](https://www.yourkit.com/images/yklogo.png)

YourKit supports open source projects with innovative and intelligent tools for monitoring and profiling Java and .NET applications.
YourKit is the creator of [YourKit Java Profiler](https://www.yourkit.com/java/profiler/), [YourKit .NET Profiler](https://www.yourkit.com/.net/profiler/), [YourKit YouMonitor](https://www.yourkit.com/youmonitor/).


