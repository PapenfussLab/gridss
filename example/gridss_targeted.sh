#!/bin/bash
#
# Runs GRIDSS on a set of targeted regions
#
# Note that to properly call compound breakpoints, this script should be run
# twice: once for the intial GRIDSS calls, then again on the GRIDSS calls
# This enables accurate calling of complex events for which only some of the
# breakpoints involved fall within the initially targeted region
#
TARGETING_INPUT=regions_to_call.bed
INPUT=chr12.1527326.DEL1024.bam
REFERENCE=hg19.fa
OUTPUT=${INPUT/.bam/.targeted.sv.vcf}
ASSEMBLY=${OUTPUT/.sv.vcf/.targeted.gridss.assembly.bam}
GRIDSS_JAR=../target/gridss-2.1.0-gridss-jar-with-dependencies.jar
WORKING_DIR=./
FLANKING_BASES=2000 # must be greater than 99.5% of fragments in the library

if [[ "$TARGETING_INPUT" == *.bed ]] ; then
	IN_BED=$1
elif
	# TODO VCF parsing
	echo "$TARGETING_INPUT is not a bed file" 2>&1
	exit 1
fi

JVM_ARGS="-ea \
	-Dreference_fasta="$REFERENCE" \
	-Dsamjdk.create_index=true \
	-Dsamjdk.use_async_io_read_samtools=true \
	-Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dgridss.gridss.output_to_temp_file=true \
	-cp $GRIDSS_JAR "

IN_BED_SLOPPED=$IN_BED.slopped.bed

bedtools slop -g ${REFERENCE}.fai -b $FLANKING_BASES -i - < $IN_BED
bedtools sort -i - | \
bedtools merge -d $((2 * FLANKING_BASES)) -i - > $IN_BED

# for each input file
java 
