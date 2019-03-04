#!/bin/bash
#
# Runs GRIDSS on a set of targeted regions
#
#
# WARNING: reference read counts will only be correct in targeted regions.
# If a breakpoint goes from a targeted region to a region not in the BED of
# interest then the untargeted breakend will be reported as homozygous variant
# even when it is not (since the reference-supporting reads were not included).
# This can be fixed by rerunning gridss.AnnotateReferenceReads against the
# original input file directly, or as explained below.
#
# To properly call compound breakpoints, this script should be run
# twice: once for the intial GRIDSS calls, then again on the GRIDSS calls
# This enables accurate calling of complex events for which only some of the
# breakpoints involved fall within the initially targeted region, as well
# as fixing the above issue with reference read counts.
#

BED_REGIONS_OF_INTEREST=example_regions_of_interest.bed
# If you have a VCF file, you can use the following script convert to a BED
# file that correctly handle symbolic SVs
#VCF_REGIONS_OF_INTEREST=targetted_regions.vcf # If your input is a VCF you'll need to convert it to a BED file
#Rscript gridss_targeted_vcf_to_region_bed.R --input $VCF_REGIONS_OF_INTEREST --ouput $BED_REGIONS_OF_INTEREST

REGIONS_OF_INTEREST=$BED_REGIONS_OF_INTEREST
REFERENCE=hg19.fa
INPUT=chr12.1527326.DEL1024.bam
OUTPUT=${INPUT/.bam/.targeted.sv.vcf}
ASSEMBLY=${INPUT/.bam/.targeted.assembly.bam}
GRIDSS_JAR=../target/gridss-2.2.1-gridss-jar-with-dependencies.jar

# Extract the reads flanking each breakpoint, not just overlapping.
# 2k is chosen since it it (should) be larger than the fragment size
# distribution and provides enough reference-supporting read pairs that
# bias in insert size distribution of the extract subset isn't too different
# to that of the full file.
REGION_PADDING_SIZE=2000

# Number of reads for which to calculate metrics for
STOP_METRICS_AFTER=10000000
JVM_ARGS="-ea \
	-Dreference_fasta="$REFERENCE" \
	-Dsamjdk.create_index=true \
	-Dsamjdk.use_async_io_read_samtools=true \
	-Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dsamjdk.buffer_size=4194304 \
	-Dgridss.gridss.output_to_temp_file=true \
	-cp $GRIDSS_JAR "
# for each input file
for F in $INPUT ; do
	TARGETED_BAM=$(dirname $INPUT)/$(basename $INPUT .bam).targeted.bam
	GRIDSS_WORKING_DIR=./$(basename $TARGETED_BAM).gridss.working
	TMP_DIR=./$(basename $TARGETED_BAM).gridss.working
	echo "Calculate metrics for $F based on first $STOP_METRICS_AFTER reads"
	# estimate metrics
	mkdir -p $GRIDSS_WORKING_DIR
	java -Xmx8g $JVM_ARGS gridss.analysis.CollectGridssMetrics \
			TMP_DIR=$GRIDSS_WORKING_DIR \
			ASSUME_SORTED=true \
			I=$F \
			O=$GRIDSS_WORKING_DIR/$(basename $TARGETED_BAM) \
			THRESHOLD_COVERAGE=25000 \
			FILE_EXTENSION=null \
			GRIDSS_PROGRAM=null \
			GRIDSS_PROGRAM=CollectCigarMetrics \
			GRIDSS_PROGRAM=CollectMapqMetrics \
			GRIDSS_PROGRAM=CollectTagMetrics \
			GRIDSS_PROGRAM=CollectIdsvMetrics \
			GRIDSS_PROGRAM=ReportThresholdCoverage \
			PROGRAM=null \
			PROGRAM=CollectInsertSizeMetrics \
			STOP_AFTER=$STOP_METRICS_AFTER \
			| tee log.$F.collectgridssmetrics.log 
	java -Xmx8g $JVM_ARGS gridss.analysis.CollectStructuralVariantReadMetrics \
			TMP_DIR=$TMP_DIR \
			I=$F \
			OUTPUT=$GRIDSS_WORKING_DIR/$(basename $TARGETED_BAM).sv_metrics \
			INSERT_SIZE_METRICS=$GRIDSS_WORKING_DIR/$(basename $TARGETED_BAM).insert_size_metrics \
			STOP_AFTER=$STOP_METRICS_AFTER \
			| tee log.$F.collectsvmetrics.log 
	echo "Extracting all reads from fragments overlapping $REGIONS_OF_INTEREST from $F"
	
	# IndexedExtractFullReads is much faster than ExtractFullReads when
	# using a small region of interest. Unfortunately, if the regions of
	# interest are all SVs, the insert size distribution will be biased
	## Don't perform async read-ahead for gridss.IndexedExtractFullReads
	## as there's not point in performing read-ahead of reads we're going
	## ignore anyway
	JVM_ARGS_SINGLE_THREADED="-ea \
			-Dreference_fasta="$REFERENCE" \
			-Dsamjdk.create_index=true \
			-Dsamjdk.use_async_io_write_samtools=true \
			-Dgridss.gridss.output_to_temp_file=true \
			-cp $GRIDSS_JAR "
	java -Xmx8g $JVM_ARGS_SINGLE_THREADED \
		gridss.IndexedExtractFullReads \
		B=$REGIONS_OF_INTEREST \
		REGION_PADDING_SIZE=$REGION_PADDING_SIZE \
		I=$F \
		O=$TARGETED_BAM 2>&1 | tee log.$F.extract.log &
	
	# If you are extracting a large portion of genome then
	# you should perform a full pass of the input file as it will
	# be faster than doing millions of index seeks()
	#java -Xmx8g $JVM_ARGS gridss.CollectGridssMetricsAndExtractFullReads \
	#	REGION_PADDING_SIZE=$REGION_PADDING_SIZE \
	#	REGION_BED=$REGIONS_OF_INTEREST \
	#	I=$F \
	#	O=$TARGETED_BAM 2>&1 | tee log.$F.extract.log &
done
wait

# Run GRIDSS on the extracted files
# warning: if you haven't increased your OS file ulimit above 4096 then will probably crash. Refer to the GRIDSS README for more details
java -Xmx31g $JVM_ARGS \
	-cp $GRIDSS_JAR gridss.CallVariants \
	TMP_DIR=. \
	WORKING_DIR=. \
	REFERENCE_SEQUENCE=$REFERENCE \
	INPUT=$(dirname $INPUT)/$(basename $INPUT .bam).targeted.bam \
	OUTPUT=$OUTPUT \
	ASSEMBLY=$ASSEMBLY \
	WORKER_THREADS=16 \
	2>&1 | tee -a log.gridss.$HOSTNAME.$$.log
