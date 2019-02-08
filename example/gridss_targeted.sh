#!/bin/bash
#
# Runs GRIDSS on a set of targeted regions
#
# Note that to properly call compound breakpoints, this script should be run
# twice: once for the intial GRIDSS calls, then again on the GRIDSS calls
# This enables accurate calling of complex events for which only some of the
# breakpoints involved fall within the initially targeted region
#

BED_REGIONS_OF_INTEREST=targetted_regions.bed # If it's already a bed then you don't need the next two lines
VCF_REGIONS_OF_INTEREST=targetted_regions.vcf # If your input is a VCF you'll need to convert it to a BED file
Rscript gridss_targeted_vcf_to_region_bed.R --input $VCF_REGIONS_OF_INTEREST --ouput #BED_REGIONS_OF_INTEREST

REGIONS_OF_INTEREST=colo829.region.bed
REFERENCE=/data/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta
OUTPUT=colo829.lite.vcf
ASSEMBLY=assembly.lite.bam
GRIDSS_JAR=gridss-2.1.2-gridss-jar-with-dependencies.jar

# This example uses a paired tumour/normal
NORMAL=COLO829R_dedup.realigned.bam
TUMOUR=COLO829T_dedup.realigned.bam

# Extract the reads flanking each breakpoint, not just overlapping.
# 2k is chosen since it it (should) be larger than the fragment size
# distribution and provides enough reference-supporting read pairs that
# bias in insert size distribution of the extract subset isn't too different
# to that of the full file.
REGION_PADDING_SIZE=2000





INPUT=../../temp/test.bam
REFERENCE=~/hartwig/temp/Homo_sapiens.GRCh37.GATK.illumina.fa
OUTPUT=${INPUT/.bam/.targeted.sv.vcf}
ASSEMBLY=${OUTPUT/.sv.vcf/.targeted.gridss.assembly.bam}
GRIDSS_JAR=../target/gridss-2.1.2-gridss-jar-with-dependencies.jar
WORKING_DIR=../../temp/
TMP_DIR=../../temp/

JVM_ARGS="-ea \
	-Dreference_fasta="$REFERENCE" \
	-Dsamjdk.create_index=true \
	-Dsamjdk.use_async_io_read_samtools=true \
	-Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dsamjdk.buffer_size=4194304 \
	-Dgridss.gridss.output_to_temp_file=true \
	-cp $GRIDSS_JAR "
# Don't perform async read-ahead for gridss.IndexedExtractFullReads
# as there's not point in performing read-ahead of reads we're going
# ignore anyway
JVM_ARGS_SINGLE_THREADED="-ea \
        -Dreference_fasta="$REFERENCE" \
        -Dsamjdk.create_index=true \
        -Dsamjdk.use_async_io_write_samtools=true \
        -Dsamjdk.use_async_io_write_tribble=true \
        -Dgridss.gridss.output_to_temp_file=true \
        -cp $GRIDSS_JAR "


for F in $NORMAL $TUMOUR ; do
	java -Xmx8g $JVM_ARGS_SINGLE_THREADED gridss.IndexedExtractFullReads B=$REGIONS_OF_INTEREST REGION_PADDING_SIZE=$REGION_PADDING_SIZE I=$F O=$F.lite.bam 2>&1 | tee log.$F.extract.log &
	# If you have are extracting a large portion of genome then
	# you should perform a full pass of the input file as it will
	# be faster than doing millions of index seeks()
	#java -Xmx8g $JVM_ARGS gridss.ExtractFullReads REGION_PADDING_SIZE=$REGION_PADDING_SIZE B=$REGIONS_OF_INTEREST I=$F O=$F.lite.bam 2>&1 | tee log.$F.extract.log &
done
wait
# Run GRIDSS on the extracted files
# warning: if you haven't increased your OS file ulimit above 4096 then will probably crash.
# 
java -Xmx31g $JVM_ARGS \
	-cp $GRIDSS_JAR gridss.CallVariants \
	TMP_DIR=. \
	WORKING_DIR=. \
	REFERENCE_SEQUENCE=$REFERENCE \
	INPUT=$NORMAL.lite.bam \
	INPUT=$TUMOUR.lite.bam \
	OUTPUT=$OUTPUT \
	ASSEMBLY=$ASSEMBLY \
	WORKER_THREADS=16 \
	2>&1 | tee -a log.gridss.$HOSTNAME.$$.log
