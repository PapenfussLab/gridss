#!/bin/bash
#
#
. common.sh
CALLER=gridss/0.10.0-SNAPSHOT
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
	for CX_ASSEMBLY_METHOD in $GRIDSS_METHODS ; do
		for CX_K in $GRIDSS_KMER ; do
			for CX_MODEL in $GRIDSS_MODEL ; do
				for CX_EXCLUSION in $GRIDSS_EXCLUSION ; do
					CX_CALLER_ARGS="$CX_K,$CX_ASSEMBLY_METHOD,$CX_MODEL,$CX_EXCLUSION"
					cx_save
					XC_OUTPUT=$CX.vcf
					XC_SCRIPT="
rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX;
cat > $CX/gridss_custom.properties << EOF
visualisation.assemblyProgress = true
visualisation.buffers = true
visualisation.bufferTrackingItervalInSeconds = 10
assembly.method=$CX_ASSEMBLY_METHOD
assembly.k=$CX_K
scoring.model=$CX_MODEL
EOF
if [[ \"$CX_EXCLUSION\" == \"SC\" ]] ; then
	echo scoring.exclude=SplitRead >> $CX/gridss_custom.properties
	echo scoring.exclude=Indel >> $CX/gridss_custom.properties
	echo scoring.exclude=SoftClip >> $CX/gridss_custom.properties
fi
if [[ \"$CX_EXCLUSION\" == \"RP\" ]] ; then
	echo scoring.exclude=UnmappedMate >> $CX/gridss_custom.properties
	echo scoring.exclude=DiscordantPair >> $CX/gridss_custom.properties
	echo scoring.exclude=SoftClip >> $CX/gridss_custom.properties
fi
if [[ \"$CX_MODEL\" == \"ReadCount\" ]] ; then
	echo variantcalling.minScore=2 >> $CX/gridss_custom.properties
	echo variantcalling.lowQuality=13 >> $CX/gridss_custom.properties
fi
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
	PER_CHR=$GRIDSS_PER_CHR \
	VERBOSITY=DEBUG \
	BLACKLIST=$BLACKLIST \
	&& \
cp $CX/breakend.vcf $CX.vcf
		"
					xc_exec
				done
			done
		done
	done
done

