#!/bin/bash
#
# Example gridss pipeline
#
NORMAL=normal.bam
TUMOUR=tumour.bam
BLACKLIST=wgEncodeDacMapabilityConsensusExcludable.bed
REFERENCE=~/reference_genomes/human/hg19.fa
OUTPUT=somatic.gridss.vcf
GRIDSS_JAR=~/target/gridss-0.11.5-jar-with-dependencies.jar

if [[ ! -f "$INPUT" ]] ; then
	echo "Missing $INPUT input file."
	exit 1
fi
if ! which bowtie2 >/dev/null 2>&1 ; then
	echo "Missing bowtie2. Please add to PATH"
	exit 1
fi
if [[ ! -f "$REFERENCE" ]] ; then
	echo "Missing reference genome $REFERENCE. Update the REFERENCE variable in the shell script to your hg19 location"
	exit 1
fi
if [[ ! -f "$REFERENCE.1.bt2" ]] ; then
	echo "Missing bowtie2 index for $REFERENCE. Could not find $REFERENCE.1.bt2. Create an index or symlink the 6 index files to the expected file names."
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
	INPUT="$NORMAL" IC=1 \
	INPUT="$TUMOUR" IC=2 \
	OUTPUT="$OUTPUT" \
	BLACKLIST="$BLACKLIST" \
	2>&1 | tee -a gridss.$HOSTNAME.$$.log

Rscript somatic.R
