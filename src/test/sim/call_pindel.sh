#!/bin/bash
#
# runs breakdancer against bams
#
. common.sh
CALLER=pindel/0.2.5b6
for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=$CALLER
	cx_save
	XC_MULTICORE=1
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="module add $CALLER ; rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	pindel -T $XC_CORES -f $CX_REFERENCE -i <( echo $CX_BAM $CX_READ_FRAGMENT_LENGTH sample ) -o $CX/out.pindel && \
	pindel2vcf -co 2 -P $CX/out.pindel -r $CX_REFERENCE -R hs37d5 -d 20130330 -v $XC_OUTPUT
	"
	xc_exec
done

