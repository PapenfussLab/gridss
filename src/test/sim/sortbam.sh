#!/bin/bash
#
# Sorts bam files
#
. common.sh

for BAM in $DATA_DIR/*.su.bam ; do
	cx_load $BAM
	# So we don't overwrite the output from the alignment itself
	XC_SUFFIX=.sort
	XC_OUTPUT=$CX.sc.bam

	SORT_PICARD=" SortSam I=$CX.sq.tmp.bam O=$CX.sc.tmp.bam SO=coordinate && \ 
	samtools index $CX.sc.tmp.bam "
	SORT_NOVO=" novosort -c $XC_CORES -i -o $CX.sc.tmp.bam $CX.sq.tmp.bam "
	XC_SCRIPT="
	FixMateInformation I=$CX.su.bam O=$CX.sq.tmp.bam SORT_ORDER=queryname
	ValidateSamFile I=$CX.su.tmp.bam > $CX.sam.validation 2>&1
	$SORT_NOVO && \
	mv $CX.sc.tmp.bam.bai $CX.sc.bam.bai && \
	mv $CX.sq.tmp.bam $CX.sq.bam && \
	mv $CX.sc.tmp.bam $CX.sc.bam
	"
	xc_exec
done
