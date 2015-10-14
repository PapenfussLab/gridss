#!/bin/bash
#
# runs variantionhunter
#
. common.sh

CL_DIR=$BASE_DIR/tools/variationhunter/CommonLawRelease
HG19_DIR=$CL_DIR/Hg19_NecessaryFiles
PATH=$CL_DIR/selection:$CL_DIR/clustering:$PATH

for DIVET in $DATA_DIR/*.output_DIVET.vh ; do
	cx_load $DIVET
	CX_DIVET=$DIVET
	CX_CALLER=variationhunter
	cx_save
	XC_OUTPUT=$CX.out.SV
	cat > $CX.lib << EOF
1
lib sample $DIVET $CX_ALIGNER_MIN_CONCORD $CX_ALIGNER_MAX_CONCORD $CX_READ_LENGTH
EOF
	# -s and -f parameters do nothing but are still required args
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	VH \
		--chro $HG19_DIR/AllChro \
		--init $HG19_DIR/initInfo \
		--lib $CX.lib \
		--repeat $HG19_DIR/Hg19.Satellite \
		--gap $HG19_DIR/hg19_Gap.Table.USCS.Clean \
		--output $CX.cluster.out \
		--outputRead $CX.name \
		--maxmapping 500 \
		--prunprob 0.001 && \
	setCover \
		-l $CX.lib \
		-r $CX.name \
		-c $CX.cluster.out \
		-t 10000 \
		-o $CX.out.SV
	# TODO VH -> VCF
	"
	xc_exec
done

