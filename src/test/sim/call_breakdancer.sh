#!/bin/bash
#
# runs breakdancer against bams
#
. common.sh
CALLER=breakdancer/1.4.5
for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=$CALLER
	CX_BAM2CFG_FLAGS=-m
	CX_CALLER_FLAGS=""
	cx_save
	XC_OUTPUT=$CX.vcf
	echo "readgroup:1	platform:illumina	map:$CX_BAM	readlen:$CX_READ_LENGTH	mean:$CX_READ_FRAGMENT_LENGTH	std:$CX_READ_FRAGMENT_STDDEV	exe:samtools view" > $CX.cfg
	XC_SCRIPT="module add samtools $CALLER ; rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	ln -s $BAM input.bam
	bam2cfg.pl $CX_BAM2CFG_FLAGS $BAM > input.cfg
	breakdancer-max $CX_CALLER_FLAGS input.cfg > ouput.ctx
	$BASE_DIR/breakdancer2vcf.py $CX_READ_FRAGMENT_LENGTH < ouput.ctx > $XC_OUTPUT
	"
	xc_exec
done

