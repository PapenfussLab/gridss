#!/bin/bash
#
# runs breakdancer against bams
#
. common.sh
PATH=$PATH

for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=gatk
	cx_save
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	gatk -T HaplotypeCaller -R $CX_REFERENCE -I $BAM -o $XC_OUTPUT
	"
	xc_exec
done

