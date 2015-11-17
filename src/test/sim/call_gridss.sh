#!/bin/bash
#
#
. common.sh
CALLER=gridss/0.9.0
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
	BLACKLIST=$CX_BLACKLIST
	if [[ "$BLACKLIST" == "" ]] ; then
		BLACKLIST=null
	fi
	CX_BAM=$BAM
	CX_CALLER=$CALLER
	#CX_CALLER_ARGS="assembly.method=Subgraph"
	cx_save
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX;
cat > $CX/gridss_custom.properties << EOF
#variantcalling.minIndelSize = 1
visualisation.assemblyProgress = true
visualisation.buffers = true
visualisation.bufferTrackingItervalInSeconds = 10
$CX_CALLER_ARGS
EOF
	java \
		-ea \
		-Xmx64g \
		-cp $FULL_JAR \
		au.edu.wehi.idsv.Idsv \
		TMP_DIR=$CX \
		WORKING_DIR=$CX \
		INPUT=$BAM \
		INPUT_MIN_FRAGMENT_SIZE=$(( CX_READ_FRAGMENT_LENGTH - 3 * CX_READ_FRAGMENT_STDDEV )) \
		INPUT_MAX_FRAGMENT_SIZE=$(( CX_READ_FRAGMENT_LENGTH + 3 * CX_READ_FRAGMENT_STDDEV )) \
		OUTPUT=$CX/breakend.vcf \
		REFERENCE=$CX_REFERENCE \
		CONFIGURATION_FILE=$CX/gridss_custom.properties \
		PER_CHR=false \
		VERBOSITY=DEBUG \
		BLACKLIST=$BLACKLIST \
		&& \
	cp $CX/breakend.vcf $CX.vcf



	"
	xc_exec
done

