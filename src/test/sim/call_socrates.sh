#!/bin/bash
#
# runs socrates against bams
#
. common.sh

for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	if [[ -f $BAM.bt2.bam ]] ; then
		echo "Using $BAM.bt2.bam as proxy"
		CX_BAM=$BAM.bt2.bam
	fi
	CX_CALLER=socrates/1.12
	CX_CALLER_ARGS=
	if [[ "$CX_ALIGNER_SOFTCLIP" == 0 ]] ; then
		echo "Socrates: skipping end-to-end aligned $BAM"
		continue
	fi
	cx_save
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	ln -s $CX_BAM $CX/input.bam
	ln -s $CX_BAM.bai $CX/input.bam.bai
	java -Xmx128g -jar ~/src/socrates/1.12/target/socrates-1.12-jar-with-dependencies.jar -t $XC_CORES $CX_REFERENCE $CX/input.bam && \
	$BASE_DIR/socrates2vcf.py < $CX/results_Socrates_paired_*.txt > $CX.vcf
	"
	xc_exec
done

