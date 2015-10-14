#!/bin/bash
#
# runs socrates against bams
#
. common.sh

for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=socrates/1.12
	CX_CALLER_ARGS=
	if [[ "$CX_ALIGNER_SOFTCLIP" == 0 ]] ; then
		echo "Socrates: skipping end-to-end aligned $BAM"
		continue
	fi
	cx_save
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	rm -f $CX.bam $CX.bam.bai 2>/dev/null
	ln -s $BAM $CX.bam
	ln -s $BAM.bai $CX.bam.bai
	java -Xmx128g -jar ~/src/socrates/1.12/target/socrates-1.12-jar-with-dependencies.jar -t $XC_CORES $CX_REFERENCE $CX.bam && \
	$BASE_DIR/socrates2vcf.py < $CX/results_Socrates_paired_*.txt > $CX.vcf
	"
	xc_exec
done

