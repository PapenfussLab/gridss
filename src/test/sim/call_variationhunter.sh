#!/bin/bash
#
# runs variantionhunter
#
. common.sh

CALLER=variationhunter/0.04
VH_DIR=$BASE_DIR/tools/variationhunter/CommonLawRelease
HG19_DIR=$VH_DIR/Hg19_NecessaryFiles
export PATH=$VH_DIR/selection:$VH_DIR/clustering:$PATH
export PATH=$BASE_DIR/tools/variationhunter/mrfast-2.6.0.1:$PATH

for FQ1 in $DATA_DIR/*.1.generated.fq ; do
	cx_load $FQ1
	CX_DIVET=$DIVET
	CX_CALLER=$CALLER
	CX_ALIGNER=mrfast/2.6.0.1
	cx_save
	XC_OUTPUT=$CX.vcf
	# -s and -f parameters do nothing but are still required args
	XC_SCRIPT="mkdir -p $CX $CX/align $CX 2>/dev/null;
if [[ ! -s $CX/align/output_DIVET.vh ]] ; then
	# run mrfast
	cd $CX/align
	mrfast --search $CX_REFERENCE --pe --seq1 $FQ1 --seq2 ${FQ1/.1./.2.} --discordant-vh --min $(( CX_READ_FRAGMENT_LENGTH - 3 * CX_READ_FRAGMENT_STDDEV )) --max $(( CX_READ_FRAGMENT_LENGTH + 3 * CX_READ_FRAGMENT_STDDEV )) || rm 
fi

# run VariationHunter
cd $CX
echo 1 > $CX/sample.lib
echo lib sample $CX/filtered.DIVET.vh $(( CX_READ_FRAGMENT_LENGTH - 3 * CX_READ_FRAGMENT_STDDEV )) $(( CX_READ_FRAGMENT_LENGTH + 3 * CX_READ_FRAGMENT_STDDEV )) $CX_READ_LENGTH >> $CX/sample.lib

awk -f $BASE_DIR/variationhunter_divit_filter.awk $CX/align/output_DIVET.vh > $CX/filtered.DIVET.vh

VH \
	--chro $HG19_DIR/AllChro \
	--init $HG19_DIR/initInfo \
	--lib $CX/sample.lib \
	--repeat $HG19_DIR/Hg19.Satellite \
	--gap $HG19_DIR/hg19_Gap.Table.USCS.Clean \
	--output $CX/sample.cluster.out \
	--outputRead $CX/sample.name \
	--maxmapping 500 \
	--prunprob 0.001 && \
multiInd_SetCover \
	-l $CX/sample.lib \
	-r $CX/sample.name \
	-c $CX/sample.cluster.out \
	-t 50000 \
	-o $CX/out.sv && \
$BASE_DIR/variationhunter2vcf.py < out.sv > $XC_OUTPUT
	"
	xc_exec
done

