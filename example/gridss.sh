#!/bin/bash
#
# Example gridss pipeline for single sample analysis
#
INPUT=chr12.1527326.DEL1024.bam
BLACKLIST=wgEncodeDacMapabilityConsensusExcludable.bed
REFERENCE=hg19.fa
OUTPUT=${INPUT/.bam/.sv.vcf}
ASSEMBLY=${OUTPUT/.sv.vcf/.gridss.assembly.bam}
GRIDSS_JAR=../target/gridss-2.2.3-gridss-jar-with-dependencies.jar
WORKER_THREADS=16

if [[ ! -f "$INPUT" ]] ; then
	echo "Missing $INPUT input file."
	exit 1
fi
if ! which bwa >/dev/null 2>&1 ; then
	echo "Missing bwa. Please add to PATH"
	exit 1
fi
if [[ ! -f "$REFERENCE" ]] ; then
	echo "Missing reference genome $REFERENCE. Update the REFERENCE variable in the shell script to your hg19 location"
	echo "For the example file chr12.1527326.DEL1024.bam, ReorderSam can be used to match to your version of hg19. In the case of this example, only \"chr12\" is required to exist, and difference in alternate contigs can be ignored (using ALLOW_INCOMPLETE_DICT_CONCORDANCE=true)."
	echo "For real data, please ensure that all BAM files are aligned to the same reference, and the reference supplied to GRIDSS matches that used for alignment."
	exit 1
fi
if [[ ! -f "$REFERENCE.bwt" ]] ; then
	echo "Missing bwa index for $REFERENCE. Could not find $REFERENCE.bwt. Create a bwa index (using \"bwa index $REFERENCE\") or symlink the index files to the expected file names."
	exit 1
fi
if [[ ! -f $GRIDSS_JAR ]] ; then
	echo "Missing $GRIDSS_JAR. Update the GRIDSS_JAR variable in the shell script to your location"
	exit 1
fi
if ! which java >/dev/null 2>&1 ; then
	echo "Missing java. Please add java 1.8 or later to PATH"
	exit 1
fi
JAVA_VERSION="$(java -version 2>&1 | grep version )"
if [[ ! "$JAVA_VERSION" =~ "\"1.8" ]] ; then
	echo "Detected $JAVA_VERSION. GRIDSS requires Java 1.8 or later."
	exit 1
fi
if ! which Rscript >/dev/null 2>&1 ; then
	echo "Missing R installation. Please add Rscript to PATH"
	exit 1
fi

ulimit -n $(ulimit -Hn) # Reduce likelihood of running out of open file handles 
unset DISPLAY # Prevents errors attempting to connecting to an X server when starting the R plotting device

# -Dreference_fasta is only required for CRAM input files
# -Dgridss.gridss.output_to_temp_file=true allows GRIDSS to continue where it left off without data errors due to truncated files
# -Dsamjdk.create_index=true is required for multi-threaded operation
# -Dsamjdk.use_async_io allow for async read/write on background threads which improves BAM I/O performancce
java -ea -Xmx31g \
	-Dreference_fasta="$REFERENCE" \
	-Dsamjdk.create_index=true \
	-Dsamjdk.use_async_io_read_samtools=true \
	-Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dsamjdk.buffer_size=4194304 \
	-Dgridss.gridss.output_to_temp_file=true \
	-cp $GRIDSS_JAR gridss.CallVariants \
	TMP_DIR=. \
	WORKING_DIR=. \
	REFERENCE_SEQUENCE="$REFERENCE" \
	INPUT="$INPUT" \
	OUTPUT="$OUTPUT" \
	ASSEMBLY="$ASSEMBLY" \
	BLACKLIST="$BLACKLIST" \
	WORKER_THREADS=$WORKER_THREADS \
	2>&1 | tee -a gridss.$HOSTNAME.$$.log

