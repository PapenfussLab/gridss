#!/bin/bash
#
# runs breakdancer against bams
#
. common.sh
PATH=$PATH

for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=samtools
	cx_save
	XC_MEMORY=1024
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	samtools mpileup -uf $CX_REFERENCE $BAM | bcftools view -v indels,other - > $XC_OUTPUT"
	xc_exec
done

