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

# If you have multiple input files, repeat for each input
java -Xmx8g $JVM_ARGS gridss.ExtractFullReads \
	B=$REGIONS_OF_INTEREST \
	I=$INPUT \
	O=$INPUT.targeted.bam \
	TMP_DIR=../../temp
	
# Then run GRIDSS on all the input files together by specifying INPUT= multiple times
java -Xmx16g $JVM_ARGS \
	-cp $GRIDSS_JAR gridss.CallVariants \
	TMP_DIR=. \
	WORKING_DIR=. \
	REFERENCE_SEQUENCE=$REFERENCE \
	INPUT=$INPUT.targeted.bam \
	OUTPUT=$OUTPUT \
	ASSEMBLY=$ASSEMBLY \
	2>&1 | tee -a log.gridss.$HOSTNAME.$$.log
