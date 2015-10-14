#!/bin/bash
#
# runs prism against bams
#
. common.sh
CALLER=varscan/2.4.0
for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=$CALLER
	CX_CALLER_FLAGS="--p-value 0.99"
	cx_save
	XC_OUTPUT=$CX.vcf
	XC_MEMORY=1024
	XC_SCRIPT="module add $CALLER ; rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	samtools mpileup -f $CX_REFERENCE $CX_BAM | \
	java -jar \\$VARSCAN_JAR mpileup2indel $CX_CALLER_FLAGS > $CX/out && \
	java -jar \\$VARSCAN_JAR fpfilter $CX/out <( bam-readcounts $CX_BAM ) --output-file $CX/pass --filtered-file $CX/fail && \
	$BASE_DIR/varscan2vcf.py < $CX/pass > $CX/pass.vcf && \
	$BASE_DIR/varscan2vcf.py < $CX/fail | sed 's/PASS/fail/g' >> $CX/fail.vcf && \
	vcf-concat $CX/pass.vcf $CX/fail.vcf > $CX.vcf
	"
	xc_exec
	exit
done

