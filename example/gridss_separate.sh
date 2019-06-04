#!/bin/bash
#
# Example gridss pipeline for single sample analysis
#
INPUT=chr12.1527326.DEL1024.bam
BLACKLIST=wgEncodeDacMapabilityConsensusExcludable.bed
REFERENCE=hg19.fa
OUTPUT=${INPUT/.bam/.sv.vcf}
ASSEMBLY=${OUTPUT/.sv.vcf/.gridss.assembly.bam}
GRIDSS_JAR=../target/gridss-2.3.0-gridss-jar-with-dependencies.jar
WORKING_DIR=.


THREADS=$(nproc)

if [[ ! -f "$INPUT" ]] ; then
	echo "Missing $INPUT input file."
	exit 1
fi
if ! which bwa >/dev/null 2>&1 ; then
	echo "Missing bwa. Please add to PATH"
	exit 1
fi
if [[ ! -f "$REFERENCE" ]] ; then
	echo "Missing reference genome $REFERENCE. Update the REFERENCE variable in the shell script to your hg19 location"
	echo "For the example file chr12.1527326.DEL1024.bam, ReorderSam can be used to match to your version of hg19. In the case of this example, only \"chr12\" is required to exist, and difference in alternate contigs can be ignored (using ALLOW_INCOMPLETE_DICT_CONCORDANCE=true)."
	echo "For real data, please ensure that all BAM files are aligned to the same reference, and the reference supplied to GRIDSS matches that used for alignment."
	exit 1
fi
if [[ ! -f "$REFERENCE.bwt" ]] ; then
	echo "Missing bwa index for $REFERENCE. Could not find $REFERENCE.bwt. Create a bwa index (using \"bwa index $REFERENCE\") or symlink the index files to the expected file names."
	exit 1
fi
if [[ ! -f $GRIDSS_JAR ]] ; then
	echo "Missing $GRIDSS_JAR. Update the GRIDSS_JAR variable in the shell script to your location"
	exit 1
fi
if ! which java >/dev/null 2>&1 ; then
	echo "Missing java. Please add java 1.8 or later to PATH"
	exit 1
fi
JAVA_VERSION="$(java -version 2>&1 | grep version )"
if [[ ! "$JAVA_VERSION" =~ "\"1.8" ]] ; then
	echo "Detected $JAVA_VERSION. GRIDSS requires Java 1.8 or later."
	exit 1
fi
if ! which Rscript >/dev/null 2>&1 ; then
	echo "Missing R installation. Please add Rscript to PATH"
	exit 1
fi

ulimit -n $(ulimit -Hn) # Reduce likelihood of running out of open file handles 
# Prevents errors attempting to connecting to an X server when starting the R plotting device
unset DISPLAY
mkdir -p $WORKING_DIR

# -Dreference_fasta=$REFERENCE is only required for CRAM input files
# -Dgridss.gridss.output_to_temp_file=true allows GRIDSS to continue where it left off without data errors due to truncated files
# -Dsamjdk.create_index=true is required for multi-threaded operation
# -Dsamjdk.use_async_io allow for async read/write on background threads which improves BAM I/O performancce
JVM_ARGS="-ea \
	-Dsamjdk.create_index=true \
	-Dsamjdk.use_async_io_read_samtools=true \
	-Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dgridss.gridss.output_to_temp_file=true \
	-cp $GRIDSS_JAR"


################
# input file preprocessing - this must be done for each input file
#
INPUT_WORKING_DIR=$WORKING_DIR/$(basename $INPUT).gridss.working
INPUT_SV_BAM=$INPUT_WORKING_DIR/$(basename $INPUT).sv.bam
INPUT_EXTRACTED_BAM=$INPUT_WORKING_DIR/gridss.tmp.querysorted.$(basename $INPUT).sv.bam
INPUT_TAGGED_BAM=$INPUT_WORKING_DIR/gridss.tmp.withtags.$(basename $INPUT).sv.bam
if [[ ! -f $INPUT_SV_BAM ]] ; then
	mkdir -p $INPUT_WORKING_DIR
	if [[ ! -f $INPUT_WORKING_DIR/$(basename $INPUT).idsv_metrics ]] ; then
		java -Xmx256m $JVM_ARGS gridss.analysis.CollectGridssMetrics \
			ASSUME_SORTED=true \
			I=$INPUT \
			O=$INPUT_WORKING_DIR/$(basename $INPUT) \
			THRESHOLD_COVERAGE=10000 \
			FILE_EXTENSION=null \
			GRIDSS_PROGRAM=null \
			GRIDSS_PROGRAM=CollectCigarMetrics \
			GRIDSS_PROGRAM=CollectMapqMetrics \
			GRIDSS_PROGRAM=CollectTagMetrics \
			GRIDSS_PROGRAM=CollectIdsvMetrics \
			GRIDSS_PROGRAM=ReportThresholdCoverage \
			PROGRAM=null \
			PROGRAM=CollectInsertSizeMetrics \
			PROGRAM=CollectAlignmentSummaryMetrics \
			PROGRAM=QualityScoreDistribution
		if [[ $? -ne 0 ]] ; then
			echo "gridss.analysis.CollectGridssMetrics returned non-zero exit status - aborting"
			exit $?
		fi
	fi
	if [[ ! -f $INPUT_EXTRACTED_BAM ]] ; then
		java -Xmx1G $JVM_ARGS gridss.ExtractSVReads \
			REFERENCE_SEQUENCE=$REFERENCE \
			I=$INPUT \
			O=/dev/stdout \
			METRICS_OUTPUT=$INPUT_WORKING_DIR/$(basename $INPUT).sv_metrics \
			INSERT_SIZE_METRICS=$INPUT_WORKING_DIR/$(basename $INPUT).insert_size_metrics \
			UNMAPPED_READS=false \
			INCLUDE_DUPLICATES=true \
			MIN_CLIP_LENGTH=5 \
		| samtools sort -O bam -T $INPUT_WORKING_DIR/samtools.sort.n.tmp -n -l 1 -@ $THREADS -o $INPUT_EXTRACTED_BAM -
		if [[ ${PIPESTATUS[0]} -ne 0 ]] ; then
			echo "gridss.ExtractSVReads returned non-zero exit status - aborting"
			rm $INPUT_EXTRACTED_BAM
			exit ${PIPESTATUS[0]}
		fi
		if [[  ${PIPESTATUS[1]} -ne 0 ]] ; then
			echo "samtools sort returned non-zero exit status - aborting"
			rm $INPUT_EXTRACTED_BAM
			exit ${PIPESTATUS[1]}
		fi
	fi
	
	if [[ ! -f $INPUT_TAGGED_BAM ]] ; then
		java -Xmx1G $JVM_ARGS gridss.ComputeSamTags \
			TMP_DIR=$WORKING_DIR \
			WORKING_DIR=$WORKING_DIR \
			REFERENCE_SEQUENCE=$REFERENCE \
			COMPRESSION_LEVEL=0 \
			ASSUME_SORTED=true \
			I=$INPUT_EXTRACTED_BAM \
			O=/dev/stdout \
		| samtools sort -O bam -T $INPUT_WORKING_DIR/samtools.sort.c.tmp -@ $THREADS -o $INPUT_TAGGED_BAM -
		if [[ ${PIPESTATUS[0]} -ne 0 ]] ; then
			echo "gridss.ComputeSamTags returned non-zero exit status - aborting"
			rm $INPUT_TAGGED_BAM
			exit ${PIPESTATUS[0]}
		fi
		if [[  ${PIPESTATUS[1]} -ne 0 ]] ; then
			echo "samtools sort returned non-zero exit status - aborting"
			rm $INPUT_TAGGED_BAM
			exit ${PIPESTATUS[1]}
		fi
		rm $INPUT_EXTRACTED_BAM
	fi
	
	java -Xmx3G $JVM_ARGS gridss.SoftClipsToSplitReads \
		TMP_DIR=$WORKING_DIR \
		WORKING_DIR=$WORKING_DIR \
		REFERENCE_SEQUENCE=$REFERENCE \
		I=$INPUT_TAGGED_BAM \
		O=$INPUT_SV_BAM \
		WORKER_THREADS=$THREADS \
		ALIGNER_COMMAND_LINE=null \
		ALIGNER_COMMAND_LINE='bwa' \
		ALIGNER_COMMAND_LINE='mem' \
		ALIGNER_COMMAND_LINE='-t' \
		ALIGNER_COMMAND_LINE='%3$d' \
		ALIGNER_COMMAND_LINE='%2$s' \
		ALIGNER_COMMAND_LINE='%1$s'
	# or if you prefer to use bowtie2
	#ALIGNER_COMMAND_LINE='bowtie2' \
	#ALIGNER_COMMAND_LINE='--threads' \
	#ALIGNER_COMMAND_LINE='%3$d' \
	#ALIGNER_COMMAND_LINE='--local' \
	#ALIGNER_COMMAND_LINE='--mm' \
	#ALIGNER_COMMAND_LINE='--reorder' \
	#ALIGNER_COMMAND_LINE='-x' \
	#ALIGNER_COMMAND_LINE='%2$s' \
	#ALIGNER_COMMAND_LINE='-U' \
	#ALIGNER_COMMAND_LINE='%1$s'
	if [[ $? -ne 0 ]] ; then
		echo "gridss.SoftClipsToSplitReads returned non-zero exit status - aborting"
		rm $INPUT_SV_BAM
		exit $?
	fi
	rm $INPUT_TAGGED_BAM
fi


################
# Assembly
#
if [[ ! -f $ASSEMBLY ]] ; then
	java -Xmx31G $JVM_ARGS gridss.AssembleBreakends \
		TMP_DIR=$WORKING_DIR \
		WORKING_DIR=$WORKING_DIR \
		REFERENCE_SEQUENCE=$REFERENCE \
		INPUT=$INPUT \
		OUTPUT=$ASSEMBLY \
		WORKER_THREADS=$THREADS \
		BLACKLIST=$BLACKLIST
	if [[ $? -ne 0 ]] ; then
		echo "gridss.AssembleBreakends returned non-zero exit status - aborting"
		rm $ASSEMBLY
		exit $?
	fi
fi
if [[ ! -f $WORKING_DIR/$(basename $ASSEMBLY).gridss.working/$(basename $ASSEMBLY).idsv_metrics ]] ; then
	java -Xmx256m $JVM_ARGS gridss.analysis.CollectGridssMetrics \
		ASSUME_SORTED=true \
		I=$ASSEMBLY \
		O=$WORKING_DIR/$(basename $ASSEMBLY).gridss.working/$(basename $ASSEMBLY) \
		THRESHOLD_COVERAGE=10000 \
		FILE_EXTENSION=null \
		GRIDSS_PROGRAM=null \
		GRIDSS_PROGRAM=CollectCigarMetrics \
		GRIDSS_PROGRAM=CollectMapqMetrics \
		GRIDSS_PROGRAM=CollectTagMetrics \
		GRIDSS_PROGRAM=CollectIdsvMetrics \
		GRIDSS_PROGRAM=ReportThresholdCoverage \
		PROGRAM=null \
		PROGRAM=CollectInsertSizeMetrics
	if [[ $? -ne 0 ]] ; then
		echo "gridss.analysis.CollectGridssMetrics returned non-zero exit status - aborting"
		exit $?
	fi
fi
if [[ ! -f $WORKING_DIR/$(basename $ASSEMBLY).gridss.working/$(basename $ASSEMBLY).sv.bam ]] ; then
	java -Xmx4G $JVM_ARGS -Dgridss.async.buffersize=16 gridss.SoftClipsToSplitReads \
		TMP_DIR=$WORKING_DIR \
		WORKING_DIR=$WORKING_DIR \
		REFERENCE_SEQUENCE=$REFERENCE \
		I=$ASSEMBLY \
		O=$WORKING_DIR/$(basename $ASSEMBLY).gridss.working/$(basename $ASSEMBLY).sv.bam \
		REALIGN_ENTIRE_READ=true \
		WORKER_THREADS=$THREADS
	if [[ $? -ne 0 ]] ; then
		echo "gridss.SoftClipsToSplitReads returned non-zero exit status - aborting"
		rm $WORKING_DIR/$(basename $ASSEMBLY).gridss.working/$(basename $ASSEMBLY).sv.bam
		exit $?
	fi
fi
################
# Variant calling
#
if [[ ! -f $OUTPUT ]] ; then
	if [[ ! -f $WORKING_DIR/$(basename $OUTPUT).unannotated.vcf ]] ; then
		java -Xmx8G $JVM_ARGS gridss.IdentifyVariants \
			TMP_DIR=$WORKING_DIR \
			WORKING_DIR=$WORKING_DIR \
			REFERENCE_SEQUENCE=$REFERENCE \
			INPUT=$INPUT \
			OUTPUT_VCF=$WORKING_DIR/$(basename $OUTPUT).unannotated.vcf \
			ASSEMBLY=$ASSEMBLY \
			WORKER_THREADS=$THREADS \
			BLACKLIST=$BLACKLIST
		if [[ $? -ne 0 ]] ; then
			echo "gridss.IdentifyVariants returned non-zero exit status - aborting"
			rm -r $WORKING_DIR/$(basename $OUTPUT).unannotated.vcf $WORKING_DIR/$(basename $OUTPUT).unannotated.vcf.gridss.working
			exit $?
		fi
	fi
	java -Xmx8G $JVM_ARGS gridss.AnnotateVariants \
		TMP_DIR=$WORKING_DIR \
		WORKING_DIR=$WORKING_DIR \
		REFERENCE_SEQUENCE=$REFERENCE \
		INPUT=$INPUT \
		INPUT_VCF=$WORKING_DIR/$(basename $OUTPUT).unannotated.vcf \
		OUTPUT_VCF=$OUTPUT \
		ASSEMBLY=$ASSEMBLY \
		WORKER_THREADS=$THREADS \
		BLACKLIST=$BLACKLIST
	if [[ $? -ne 0 ]] ; then
		echo "gridss.AnnotateVariants returned non-zero exit status - aborting"
		rm $WORKING_DIR/$(basename $OUTPUT).unannotated.vcf
		exit $?
	fi
	rm $WORKING_DIR/$(basename $OUTPUT).unannotated.vcf
fi
