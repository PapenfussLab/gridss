#!/bin/bash
#
# runs breakdancer against bams
#
. common.sh
PATH=$BASE_DIR/tools/soapindel/indel_detection.release:$PATH

for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=soapindel
	cx_save
	echo "$CX_BAM	$CX_READ_FRAGMENT_LENGTH	$((CX_READ_FRAGMENT_STDDEV * 2))	$CX_READ_LENGTH" > $CX.mapping.list
	XC_OUTPUT=$CX.vcf
	XC_MEMORY=1024
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	indel_detection.ibam.pl \
		$CX.mapping.list \
		$CX_REFERENCE.DB.SPLIT \
		-cpu $XC_CORES \
		-wd $CX
	# merge and convert to valid VCF
	for RESULT in \$(find -name '*.indel.vcf') ; do
		# strip invalid VCF headers
		grep -vE '^# ' \$RESULT > \$RESULT.corrected
	done
	if [ \$(find -name '*.indel.vcf.corrected' | wc -l) -gt 0 ] ; then
		vcf-concat \$(find -name '*.indel.vcf.corrected') > $CX.vcf || rm $CX.vcf
	fi
	"
	xc_exec
done

