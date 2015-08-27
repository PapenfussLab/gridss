#!/bin/bash
#
# Example gridss pipeline
#
INPUT=chr12.1527326.DEL1024.bam
OUTPUT=${INPUT/.bam/.vcf}
REFERENCE=~/reference_genomes/human/hg19.fa
JAVA_ARGS="-ea -Xmx16g -cp ../target/gridss-*-jar-with-dependencies.jar"
CORES=$(nproc 2>/dev/null || echo 1)

# ensure libsswjni.so can be loaded
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/../lib

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

exec_gridss() {
	if [[ -f $OUTPUT ]] ; then
		return
	fi
	rm -f realign.sh
	java $JAVA_ARGS au.edu.wehi.idsv.Idsv \
		TMP_DIR=. \
		WORKING_DIR=. \
		INPUT="$INPUT" \
		OUTPUT="$OUTPUT" \
		REFERENCE="$REFERENCE" \
		SCRIPT=realign.sh \
		VERBOSITY=INFO \
		WORKER_THREADS=$CORES \
		2>&1 | tee -a gridss.$1.log
	
	if [[ -f realign.sh ]] ; then
		chmod a+x realign.sh
		./realign.sh 2>&1 | tee -a gridss.$1.log
		exec_gridss
	fi
}

exec_gridss

if [[ -f $OUTPUT ]] ; then
	java $JAVA_ARGS au.edu.wehi.idsv.VcfBreakendToBedpe \
		INPUT="$OUTPUT" \
		OUTPUT="${OUTPUT/.vcf/.bedpe}" \
		REFERENCE="$REFERENCE"
		
fi