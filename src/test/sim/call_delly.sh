#!/bin/bash
#
# runs breakdancer against bams
#
. common.sh
CALLER=delly/0.6.8
for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=$CALLER
	cx_save
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="module add $CALLER; rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	delly -t DEL -o del.vcf -g $CX_REFERENCE $CX_BAM && \
	delly -t DUP -o dup.vcf -g $CX_REFERENCE $CX_BAM && \
	delly -t INV -o inv.vcf -g $CX_REFERENCE $CX_BAM && \
	delly -t TRA -o tra.vcf -g $CX_REFERENCE $CX_BAM && \
	vcfsort *.vcf > $XC_OUTPUT
	"
	xc_exec
done

