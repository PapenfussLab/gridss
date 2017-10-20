#!/bin/bash
#
# Example of GRIDSS integrated into a NGS pipeline
#
# This approach is more computationally efficient than
# running GRIDSS as a stand-alone program as:
# - many steps can be integrated into earlier pipeline stages
# - some steps can be skipped since we know what aligner and settings were used
#

FQ1=example.R1.fq
FQ2=example.R2.fq
REFERENCE=~/reference_genomes/human/hg19.fa
INPUT=example.input.bam
OUTPUT=example.vcf
GRIDSS_JAR=~/bin/gridss-1.4.5-SNAPSHOT-jar-with-dependencies.jar
JVM_ARGS="
	-Dsamjdk.use_async_io_read_samtools=true 
	-Dsamjdk.use_async_io_write_samtools=true 
	-Dsamjdk.use_async_io_write_tribble=true"

	
####################
# Environment/external tools sanity checks

if [[ ! -f "$PICARD" ]] ; then
	echo "The PICARD environment variable does not point to the picard tools jar. Please add \$PICARD and point it to your picard tools location"
	exit 1
fi
if ! which bwa >/dev/null 2>&1 ; then
	echo "Missing bwa. Please add to PATH"
	exit 1
fi
if [[ ! -f "$REFERENCE" ]] ; then
	echo "Missing reference genome $REFERENCE. Update the REFERENCE variable in the shell script to your hg19 location"
	echo "For the example file chr12.1527326.DEL1024.bam, ReorderSam can be used to match to your version of hg19. In the case of this example, only \"chr12\" is required to exist, and difference in alternate contigs can be ignored (using ALLOW_INCOMPLETE_DICT_CONCORDANCE=true)."
	echo "For real data, please ensure that all BAM files are aligned to the same reference, and the reference supplied to GRIDSS matched that used for alignment."
	exit 1
fi
if [[ ! -f "$REFERENCE.bwt" ]] ; then
	echo "Missing bwa index for $REFERENCE. Could not find $REFERENCE.bwt. Create a bwa index (using \"bwa index $REFERENCE\") or symlink the index files to the expected file names."
	exit 1
fi


####################
# per INPUT file steps
if [[ ! -f "$FQ1" ]] ; then
	echo "Missing $FQ1"
	exit 1
fi
if [[ ! -f "$FQ2" ]] ; then
	echo "Missing $FQ2"
	exit 1
fi
INPUT_WORKING_DIR=$(dirname $INPUT)/$(basename $INPUT).gridss.working
INPUT_WORKING_PREFIX=$INPUT_WORKING_DIR/$(basename $INPUT)

mkdir -p $INPUT_WORKING_DIR
# bwa mem / samtools rmdup version
bwa mem -t $(nproc) -Y $REFERENCE $FQ1 $FQ2 | \
samtools fixmate -m - - | \
# Skip compression since we're in a pipe
# RECALCULATE_SA_SUPPLEMENTARY not needed since we aligned with -Y and without -M
# SOFTEN_HARD_CLIPS not needed since we supplied the -Y flag to bwa
# Since we're using bwa, we only need to write the mate info tags
java $JVM_ARGS -Dsamjdk.create_index=false \
	-cp $GRIDSS_JAR gridss.ComputeSamTags \
	REFERENCE_SEQUENCE=$REFERENCE \
	WORKING_DIR=$INPUT_WORKING_DIR \
	TMP_DIR=$INPUT_WORKING_DIR \
	I=/dev/stdin \
	O=/dev/stdout \
	COMPRESSION_LEVEL=0 \
	RECALCULATE_SA_SUPPLEMENTARY=false \
	SOFTEN_HARD_CLIPS=false \
	FIX_MATE_INFORMATION=false \
	TAGS=null \
	TAGS=R2 \
	TAGS=Q2 \
	TAGS=MC \
	TAGS=MQ \
	AS=true \
	| \
java $JVM_ARGS -Dsamjdk.create_index=false -cp $GRIDSS_JAR gridss.SoftClipsToSplitReads \
	REFERENCE_SEQUENCE=$REFERENCE \
	I=/dev/stdin \
	O=- \
	ALIGNER_STREAMING=true \
	WORKER_THREADS=$(nproc) \
	COMPRESSION_LEVEL=0 \
	| \
samtools view | less
exit
# THRESHOLD_COVERAGE=gridss.properties/maxCoverage
# Using non-POSIX-compatble process substitution feature of bash to fork the output stream
tee >(java $JVM_ARGS -cp $GRIDSS_JAR gridss.analysis.CollectGridssMetrics REFERENCE_SEQUENCE=$REFERENCE I=/dev/stdin O=$INPUT_WORKING_PREFIX THRESHOLD_COVERAGE=10000) | \
samtools sort -O BAM -T $INPUT_WORKING_DIR - | \
samtools markdup - $INPUT

# TODO: piped SoftClipsToSplitReads (can ignore sort order) STREAMING_ALIGNMENT=true RETAIN_SORT_ORDER=false

# TODO: pipeline using samtools fixmate & samtools markdup
	
	