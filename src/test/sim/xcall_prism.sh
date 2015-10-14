#!/bin/bash
#
# runs prism against bams
#
. common.sh
export PRISM_PATH=$BASE_DIR/tools/prism/PRISM_1_1_6
PATH=$PATH:$PRISM_PATH/toolkit:$PRISM_PATH/bin

for BAM in $DATA_DIR/*.sq.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=prism
	CX_CALLER_FLAGS=""
	cx_save
	if ls $CX_REFERENCE.DB.SPLIT/*.fa 2>&1 > /dev/null ; then
		# found reference
		echo -n
	else
		echo "Unable to file per-chromosome fasta references in $CX_REFERENCE.DB.SPLIT. Unable to run PRISM"
		exit 1
	fi
	XC_OUTPUT=$CX.vcf
	#XC_TRAP="{ rm -f $CX/*.sam }"
	XC_SCRIPT=" mkdir $CX 2>/dev/null ; cd $CX"
	for CHR_REFERENCE in $CX_REFERENCE.DB.SPLIT/*.fa ; do
		CHR=$(basename $CHR_REFERENCE .fa)
		XC_SCRIPT="$XC_SCRIPT
		samtools view $CX_BAM | gawk '{ if (\$3==\"$CHR\" || \$7==\"$CHR\") { print \$0 } }' > $CX/$CHR.sam && \
		run_PRISM.sh \
			-m $CX_READ_FRAGMENT_LENGTH \
			-e $CX_READ_FRAGMENT_STDDEV \
			-I $CX/$CHR/in \
			-O $CX/$CHR/out \
			-r $CX_REFERENCE.DB.SPLIT/$CHR.fa \
			-i $CX/$CHR.sam \
			-a BWA \
			-l $CX_READ_LENGTH \
			$CX_CALLER_FLAGS && \
		cat \$(find . -name split_all.sam_ns_rmmul_cigar_sorted_sv) | $BASE_DIR/prism2vcf.py > $CX.vcf
		"
	done
	xc_exec
done

