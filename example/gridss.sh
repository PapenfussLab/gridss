#!/bin/bash
#
# Example gridss pipeline
#
INPUT=chr12.1527326.DEL1024.bam
BLACKLIST=wgEncodeDacMapabilityConsensusExcludable.bed
REFERENCE=~/reference_genomes/human/hg19.fa
OUTPUT=${INPUT/.bam/.gridss.vcf}
GRIDSS_JAR=~/bin/gridss-0.11.1-jar-with-dependencies.jar

if [[ ! -f "$INPUT" ]] ; then
	echo "Missing $INPUT input file."
	exit 1
fi
if ! which bwa >/dev/null 2>&1 ; then
	echo "Missing bwa executable. Please add to PATH"
	exit 1
fi
if [[ ! -f "$REFERENCE" ]] ; then
	echo "Missing reference genome $REFERENCE. Update the REFERENCE variable in the shell script to your hg19 location"
	exit 1
fi
if [[ ! -f "$REFERENCE.bwt" ]] ; then
	echo "Missing bwa index for $REFERENCE. Could not find $REFERENCE.1.bt2. Create an index (bwa index $REFERENCE $REFERENCE) or symlink the index files to the expected file names."
	exit 1
fi

if [[ ! -f $GRIDSS_JAR ]] ; then
	echo "Missing $GRIDSS_JAR. Update the GRIDSS_JAR variable in the shell script to your location"
	exit 1
fi

JAVA_VERSION="$(java -version 2>&1 | head -1)"
if [[ ! "$JAVA_VERSION" =~ "\"1.8" ]] ; then
	echo "Detected $JAVA_VERSION. GRIDSS requires Java 1.8 or later."
	exit 1
fi

java -ea -Xmx16g -cp $GRIDSS_JAR au.edu.wehi.idsv.Idsv \
	TMP_DIR=. \
	WORKING_DIR=. \
	REFERENCE="$REFERENCE" \
	INPUT="$INPUT" IC=1 \
	OUTPUT="$OUTPUT" \
	BLACKLIST="$BLACKLIST" \
	2>&1 | tee -a gridss.$$.log


if [[ -f "$OUTPUT" ]] ; then
	java -ea -Xmx16g -cp $GRIDSS_JAR au.edu.wehi.idsv.VcfBreakendToBedpe \
		INPUT="$OUTPUT" \
		OUTPUT="${OUTPUT/.vcf/.bedpe}" \
		OUTPUT_FILTERED="${OUTPUT/.vcf/.filtered.bedpe}" \
		REFERENCE="$REFERENCE"
fi
