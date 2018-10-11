#!/bin/bash
#
# Example of GRIDSS integrated into a NGS pipeline
#
# This approach is more computationally efficient than
# running GRIDSS as a stand-alone program as:
# - many steps can be integrated into earlier pipeline stages
# - some steps can be skipped since we know what aligner and settings were used
#

# environment modules
module remove samtools bwa java R
module add bwa/0.7.15 samtools/1.6 java/1.8.0_131 R/3.5.1 mbuffer
#cd ~/test/

FQ1=example.R1.fq
FQ2=example.R2.fq
REFERENCE=~/reference_genomes/human/hg19.fa
INPUT=pipelined.example.input.bam
OUTPUT=pipelined.example.sv.vcf
RAW_GRIDSS_ASSEMBLY=${OUTPUT/.sv.vcf/.gridss.assembly.bam}
GRIDSS_JAR=~/bin/gridss-2.0.1-jar-with-dependencies.jar
GRIDSS_JVM_ARGS="
	-Dsamjdk.use_async_io_read_samtools=true 
	-Dsamjdk.use_async_io_write_samtools=true 
	-Dsamjdk.use_async_io_write_tribble=true"
# Enabling assertions performs GRIDSS data consistency checks.
# Removing these can slightly improve performance
GRIDSS_JVM_ARGS="$GRIDSS_JVM_ARGS -ea" 
GRIDSS_THRESHOLD_COVERAGE=10000
LOG_PREFIX=log.gridss.pipeline.$HOSTNAME.$$

echo Starting pipelined GRIDSS example
# Optional performance logging
dstat -cdmnpry -t --vm --top-cpu --noheader --output $LOG_PREFIX.timing.log &
TIMING_PID=$!

####################
# Environment/external tools sanity checks
if [[ ! -f "$REFERENCE" ]] ; then
	echo "Missing reference genome $REFERENCE"
	exit 1
fi
if [[ ! "$(samtools --version | head -1 2>/dev/null )" =~ "1.6" ]] ; then
	# samtools markdup is in samtools version 1.6
	echo "Missing samtools v1.6. Please add to PATH"
	exit 1
fi
if ! which bwa >/dev/null 2>&1 ; then
	echo "Missing bwa. Please add to PATH"
	exit 1
fi
if [[ ! -f "$REFERENCE.bwt" ]] ; then
	echo "Missing bwa index for $REFERENCE. Creating."
	bwa index $REFERENCE $REFERENCE
fi
if [[ ! -f "$REFERENCE.fai" ]] ; then
	echo "Missing reference genome index for $REFERENCE. Creating."
	samtools faidx $REFERENCE
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

### bwa alignment with samtools 1.6 duplicate marking
if [[ ! -f $INPUT ]] ; then
	if [[ ! -f $INPUT_WORKING_PREFIX.unsorted.bam ]] ; then
		# First step of pipeline: alignment
		nice bwa mem -t $(nproc) -Y $REFERENCE $FQ1 $FQ2 2> $LOG_PREFIX.1.bwa.log | \
		# buffer bwa output to improve pipe performance
		mbuffer -m 256M | \
		# COMPRESSION_LEVEL: skip compression since we're in a pipe
		# RECALCULATE_SA_SUPPLEMENTARY: not needed since we aligned with -Y and without -M
		# SOFTEN_HARD_CLIPS:  not needed since we supplied the -Y flag to bwa
		# FIX_MATE_INFORMATION: not needed since we used bwa
		java -Xmx1G $GRIDSS_JVM_ARGS -Dsamjdk.create_index=false \
			-cp $GRIDSS_JAR gridss.ComputeSamTags \
			COMPRESSION_LEVEL=0 \
			REFERENCE_SEQUENCE=$REFERENCE \
			WORKING_DIR=$INPUT_WORKING_DIR \
			TMP_DIR=$INPUT_WORKING_DIR \
			I=/dev/stdin \
			O=/dev/stdout \
			RECALCULATE_SA_SUPPLEMENTARY=false \
			SOFTEN_HARD_CLIPS=false \
			FIX_MATE_INFORMATION=false \
			TAGS=null \
			TAGS=R2 \
			TAGS=Q2 \
			TAGS=MC \
			TAGS=MQ \
			AS=true \
			2> $LOG_PREFIX.1.ComputeSamTags.log | \
		# Although bwa reports split reads, you can identify additional split reads
		# by feeding the soft clipped reads back to bwa
		# Note that ALIGNER_STREAMING=true changes the order of the
		# records but that's ok here as we're going to sort it in the next
		# step
		java -Xmx2G $GRIDSS_JVM_ARGS -Dsamjdk.create_index=false \
			-cp $GRIDSS_JAR gridss.SoftClipsToSplitReads \
			COMPRESSION_LEVEL=0 \
			REFERENCE_SEQUENCE=$REFERENCE \
			I=/dev/stdin \
			O=/dev/stdout \
			ALIGNER_STREAMING=true \
			WORKER_THREADS=$(nproc) \
			2> $LOG_PREFIX.1.SoftClipsToSplitReads.log | \
		# Using non-POSIX process substitution feature of bash to fork the output stream
		# Unfortunately, the current version of HTSJDK has relatively slow compression.
		# We work around this by deferring compression to samtools in a different process.
		# As the unsorted BAM is only a temporary file, we're using fast compresison
		# bwa processes records in batches. Without buffering the output
		# bwa will stall waiting for the downstream pipe to process the output
		# before aligning the next batch
		# NB: two additional threads have been added to prevent this step bottlenecking the pipe
		tee >(samtools fixmate -@ 4 -m - $INPUT_WORKING_PREFIX.unsorted.bam 2> $LOG_PREFIX.1.fixmate.log) | \
		# THRESHOLD_COVERAGE=gridss.properties/maxCoverage
		# REFERENCE_SEQUENCE: don't supply since our metrics don't care about order and Picard crashes
		#  if we supply out of order records with a reference sequence
		java -Xmx1G $GRIDSS_JVM_ARGS -cp $GRIDSS_JAR \
			gridss.analysis.CollectGridssMetrics \
			ASSUME_SORTED=true \
			I=/dev/stdin \
			O=$INPUT_WORKING_PREFIX \
			THRESHOLD_COVERAGE=$GRIDSS_THRESHOLD_COVERAGE \
			2> $LOG_PREFIX.1.CollectGridssMetrics.log
	fi
	# Second step of pipeline: sorting
	samtools sort -@ $(nproc) -O BAM -l 0 -T $INPUT_WORKING_DIR $INPUT_WORKING_PREFIX.unsorted.bam 2> $LOG_PREFIX.2.1.sort.log | \
	samtools markdup - - -O BAM 2> $LOG_PREFIX.2.2.markdup.log | \
	# Using non-POSIX process substitution feature of bash to fork the output stream
	# Strip out R2, Q2 tags to save space. These can be rehydrated from
	# the mate information and GRIDSS only needs them for the extracted reads
	# MC/MQ could also be stripped but they do not significantly change the file size
	tee >(samtools view -b -x R2 -x Q2 - $INPUT 2> $LOG_PREFIX.2.4.striptags.log) | \
	# MIN_CLIP_LENGTH=gridss.properties/softclip.minLength
	# This step can be skipped if the R2/Q2 tags are retained. For most data sets,
	# keeping the extracted SV reads with R2/Q2 and stripping the tags from the
	# full BAM uses less space than retain the R2/Q2 tags.
	java $GRIDSS_JVM_ARGS -Dsamjdk.create_index=true \
		-cp $GRIDSS_JAR gridss.ExtractSVReads \
		REFERENCE_SEQUENCE=$REFERENCE \
		I=/dev/stdin \
		O=$INPUT_WORKING_PREFIX.sv.bam \
		METRICS_OUTPUT=$INPUT_WORKING_PREFIX.sv_metrics \
		INSERT_SIZE_METRICS=$INPUT_WORKING_PREFIX.insert_size_metrics \
		UNMAPPED_READS=false MIN_CLIP_LENGTH=5 2> $LOG_PREFIX.2.3.ExtractSVReads.log && \
	samtools index $INPUT 2> $LOG_PREFIX.2.5.index.log && \
	# the unsorted input file can now be deleted
	rm $INPUT_WORKING_PREFIX.unsorted.bam && \
	samtools index $INPUT_WORKING_PREFIX.sv.bam 2> $LOG_PREFIX.2.6.index.log
fi
####################
# Insert other steps here

# DO NOT PERFORM INDEL REALIGNMENT
# Not only does does it significantly increase the FDR of SV calls
# but it strips hard clips and makes the SA and MC tags inconsistent.


####################
# GRIDSS SV calls
#
# Preprocess of input files has been completed, so
# we just have to perform assembly, variant detection,
# and annotation
#

# Add BLACKLIST parameter if using hg19 to speed up runtime
OUTPUT_WORKING_DIR=$(dirname $OUTPUT)/$(basename $OUTPUT).gridss.working
OUTPUT_WORKING_PREFIX=$OUTPUT_WORKING_DIR/$(basename $INPUT)
# Multiple input files require INPUT= to be specified multiple times
GRIDSS_COMMON_ARGS="
	REFERENCE_SEQUENCE=$REFERENCE
	INPUT=$INPUT
	WORKER_THREADS=$(nproc)
	"
mkdir -p $OUTPUT_WORKING_DIR
mkdir -p $RAW_GRIDSS_ASSEMBLY.gridss.working
PROCESSED_GRIDSS_ASSEMBLY=$RAW_GRIDSS_ASSEMBLY.gridss.working/$(basename $RAW_GRIDSS_ASSEMBLY).sv.bam
if [[ ! -f $PROCESSED_GRIDSS_ASSEMBLY ]] ; then	
	# GRIDSS assembly
	java $GRIDSS_JVM_ARGS -Xmx31g -Dsamjdk.create_index=true \
		-cp $GRIDSS_JAR gridss.AssembleBreakends \
		$GRIDSS_COMMON_ARGS \
		OUTPUT=$RAW_GRIDSS_ASSEMBLY \
	2> $LOG_PREFIX.3.1.AssembleBreakends.log && \
	# We now need to work out which breakpoint is supported by each breakend assembly
	# We do this by aligning the soft clipped bases of breakend assembly contigs
	# This is exactly the same process as converting soft clipped reads into split reads 
	java $GRIDSS_JVM_ARGS -Dsamjdk.create_index=false \
			-cp $GRIDSS_JAR gridss.SoftClipsToSplitReads \
			REFERENCE_SEQUENCE=$REFERENCE \
			I=$RAW_GRIDSS_ASSEMBLY \
			O=/dev/stdout \
			ALIGNER_STREAMING=true \
			WORKER_THREADS=$(nproc) \
	2> $LOG_PREFIX.3.2.SoftClipsToSplitReads.log | \
	samtools sort -@ $(nproc) -O BAM -T $RAW_GRIDSS_ASSEMBLY.gridss.working - > $PROCESSED_GRIDSS_ASSEMBLY 2> $LOG_PREFIX.3.3.sort.log && \
	samtools index $PROCESSED_GRIDSS_ASSEMBLY 2> $LOG_PREFIX.3.4.index.log && \
	java $GRIDSS_JVM_ARGS -cp $GRIDSS_JAR gridss.analysis.CollectGridssMetrics \
		REFERENCE_SEQUENCE=$REFERENCE \
		I=$RAW_GRIDSS_ASSEMBLY \
		O=$RAW_GRIDSS_ASSEMBLY.gridss.working/$(basename $RAW_GRIDSS_ASSEMBLY) \
		THRESHOLD_COVERAGE=$GRIDSS_THRESHOLD_COVERAGE \
	2> $LOG_PREFIX.3.5.CollectGridssMetrics.log 
fi
if [[ ! -f $OUTPUT_WORKING_PREFIX.breakpoint.vcf ]] ; then
	# GRIDSS variant calling
	java $GRIDSS_JVM_ARGS -Xmx4g -Dsamjdk.create_index=true \
		-cp $GRIDSS_JAR gridss.IdentifyVariants \
		$GRIDSS_COMMON_ARGS \
		ASSEMBLY=$PROCESSED_GRIDSS_ASSEMBLY \
		OUTPUT_VCF=$OUTPUT_WORKING_PREFIX.breakpoint.vcf \
	2> $LOG_PREFIX.4.1.IdentifyVariants.log 
fi
if [[ ! -f $OUTPUT ]] ; then
	# GRIDSS evidence allocation and annotation
	# AnnotateVariants combines AllocateEvidence, AnnotateReferenceCoverage and AnnotateInexactHomology
	java $GRIDSS_JVM_ARGS -Xmx4g -Dsamjdk.create_index=true \
		-cp $GRIDSS_JAR gridss.AnnotateVariants \
		$GRIDSS_COMMON_ARGS \
		ASSEMBLY=$PROCESSED_GRIDSS_ASSEMBLY \
		INPUT_VCF=$OUTPUT_WORKING_PREFIX.breakpoint.vcf \
		OUTPUT_VCF=$OUTPUT \
	2> $LOG_PREFIX.5.1.AnnotateVariants.log 
fi

# Optional performance logging
kill $TIMING_PID








