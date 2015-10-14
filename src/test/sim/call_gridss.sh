#!/bin/bash
#
#
. common.sh
CALLER=gridss/0.8.1-SNAPSHOT
FULL_JAR=~/bin/${CALLER/\//-}-jar-with-dependencies.jar

for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	if [ "$CX_READ_FRAGMENT_LENGTH" == "" ] ; then
		echo "$CX missing read fragment length information: skipping"
		continue
	fi
	if [ "$CX_ALIGNER_SOFTCLIP" == 0 ] ; then
		echo "GRIDSS: skipping end-to-end aligned $BAM"
		continue
	fi

	CX_BAM=$BAM
	CX_CALLER=$CALLER
	CX_CALLER_FLAGS=""
	cx_save
	XC_OUTPUT=$CX.vcf
	# TRUTH_VCF=$CX_REFERENCE_VCF \
	# first invokation fails since realignment not done
	# don't clean up generated files for now
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX;
exec_gridss() {
	if [[ -f $CX/breakend.vcf ]] ; then
		return
	fi
	rm -f $CX/realign.sh
	java \
		-ea \
		-Xmx32g \
		-cp $FULL_JAR \
		-Dgridss.visualisation.saveall=false \
		-Dgridss.visualisation.savetimeouts=true \
		-Dgridss.visualisation.trackAssemblyProgress=true \
		au.edu.wehi.idsv.Idsv \
		TMP_DIR=$CX \
		WORKING_DIR=$CX \
		INPUT=$BAM \
		INPUT_READ_PAIR_MAX_CONCORDANT_FRAGMENT_SIZE=$(( CX_READ_FRAGMENT_LENGTH + 3 * CX_READ_FRAGMENT_STDDEV )) \
		INPUT_READ_PAIR_MIN_CONCORDANT_FRAGMENT_SIZE=$(( CX_READ_FRAGMENT_LENGTH - 3 * CX_READ_FRAGMENT_STDDEV )) \
		MIN_INDEL_SIZE=0 \
		OUTPUT=$CX/breakend.vcf \
		REFERENCE=$CX_REFERENCE \
		SCRIPT=$CX/realign.sh \
		PER_CHR=false \
		VERBOSITY=DEBUG \
		$CX_CALLER_FLAGS
	
	if [[ -f $CX/realign.sh ]] ; then
		source $CX/realign.sh
		exec_gridss
	fi
}
exec_gridss
cp $CX/breakend.vcf $CX.vcf
	"
	xc_exec
done

