#!/bin/bash
#
# runs breakdancer against bams
#
. common.sh
PATH=$PATH:$BASE_DIR/tools/Hydra-Version-0.5.3/bin:$BASE_DIR/tools/Hydra-Version-0.5.3/scripts
# https://github.com/chapmanb/bcbb/blob/master/nextgen/scripts/hydra_to_vcf.py
HYDRA2VCF=$BASE_DIR/tools/hydra_to_vcf.py

for BAM in $DATA_DIR/*.sq.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=hydra
	cx_save
	MAD=`echo "scale=4; $CX_READ_FRAGMENT_STDDEV / 1.4826" | bc`
	MLD=`echo "scale=4; 10 * $MAD" | bc`
	MNO=`echo "scale=4; 20 * $MAD" | bc`
	# follow the hydra workflow from https://code.google.com/p/hydra-sv/wiki/TypicalWorkflow
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="
	bamToBed -tag NM -i $CX_BAM | \
	pairDiscordants.py -i stdin -m hydra > $CX.bedpe &&\
	dedupDiscordants.py -i $CX.bedpe -s 3 > $CX.deduped.bedpe &&\
	hydra -in $CX.deduped.bedpe -out $CX.breaks -mld $MLD -mno $MNO &&\
	$HYDRA2VCF $CX.breaks.final $CX_REFERENCE.2bit &&\
	mv $CX.breaks.vcf $CX.vcf
	"
	xc_exec
done

