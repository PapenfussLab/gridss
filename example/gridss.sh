#!/bin/bash
#
# Example gridss pipeline
#
INPUT=chr12.1527326.DEL1024.bam
OUTPUT=${INPUT/.bam/.gridss.vcf}
REFERENCE=~/reference_genomes/human/hg19.fa
JAVA_ARGS="-ea -Xmx16g -cp ../target/gridss-*-jar-with-dependencies.jar"

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

java $JAVA_ARGS au.edu.wehi.idsv.Idsv \
	TMP_DIR=. \
	WORKING_DIR=. \
	INPUT="$INPUT" \
	OUTPUT="$OUTPUT" \
	REFERENCE="$REFERENCE" \
	VERBOSITY=INFO \
	2>&1 | tee -a gridss.$$.log

if [[ -f "$OUTPUT" ]] ; then
	java $JAVA_ARGS au.edu.wehi.idsv.VcfBreakendToBedpe \
		INPUT="$OUTPUT" \
		OUTPUT="${OUTPUT/.vcf/.bedpe}" \
		OUTPUT_FILTERED="${OUTPUT/.vcf/.filtered.bedpe}" \
		REFERENCE="$REFERENCE"
		
fi