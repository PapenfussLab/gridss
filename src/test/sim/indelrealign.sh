#!/bin/bash
#
# Performing indel realignment
#
. common.sh

for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	if [[ "$CX_INDEL_REALIGNMENT" == "" ]] ; then
		CX_INDEL_REALIGNMENT=gatk
		CX_INDEL_REALIGNMENT_SOURCE_BAM=$BAM
		cx_save
		XC_OUTPUT=$CX.sc.bam
		XC_MEMORY=4096
		XC_SCRIPT="
		gatk -T RealignerTargetCreator \
			-R $CX_REFERENCE \
			-I $BAM \
			-o $CX.intervals && \
		gatk -T IndelRealigner \
			-R $CX_REFERENCE \
			-I $BAM \
			-targetIntervals $CX.intervals \
			-o $CX.tmp.bam && \
		SortSam I=$CX.tmp.bam O=$CX.sq.tmp.bam SO=queryname && \
		samtools index $CX.tmp.bam && \
		mv $CX.tmp.bai $CX.sc.bam.bai && \
		mv $CX.sq.tmp.bam $CX.sq.bam && \
		mv $CX.tmp.bam $CX.sc.bam
		"
		xc_exec
	fi
done
