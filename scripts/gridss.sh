#!/bin/bash
#
# GRIDSS: a sensitive structural variant calling toolkit
#
# Example ../scripts/gridss.sh  -t 4 -b wgEncodeDacMapabilityConsensusExcludable.bed -r ../hg19.fa -w out -o out/gridss.full.chr12.1527326.DEL1024.vcf -a out/gridss.full.chr12.1527326.DEL1024.assembly.bam -j ../target/gridss-2.9.0-gridss-jar-with-dependencies.jar --jvmheap 8g chr12.1527326.DEL1024.bam

getopt --test
if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
	echo 'WARNING: "getopt --test"` failed in this environment.' 1>&2
	echo "WARNING: The version of getopt(1) installed on this system might not be compatible with the GRIDSS driver script." 1>&2
fi
set -o errexit -o pipefail -o noclobber -o nounset
last_command=""
current_command=""
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
trap 'echo "\"${last_command}\" command completed with exit code $?."' EXIT
#253 forcing C locale for everything
export LC_ALL=C

USAGE_MESSAGE="
Usage: gridss.sh --reference <reference.fa> --output <output.vcf.gz> --assembly <assembly.bam> [--threads n] [--jar gridss.jar] [--workingdir <directory>] [--jvmheap 25g] [--blacklist <blacklist.bed>] [--steps All|PreProcess|Assemble|Call] [--configuration gridss.properties] [--maxcoverage 50000] [--labels input1,input2,...] input1.bam [input2.bam [...]]

	-r/--reference: reference genome to use. Must have a .fai index file and a bwa index.
	-o/--output: output VCF.
	-a/--assembly: location of the GRIDSS assembly BAM. This file will be created by GRIDSS.
	-t/--threads: number of threads to use. Defaults to the number of cores available.
	-j/--jar: location of GRIDSS jar
	-w/--workingdir: directory to place GRIDSS intermediate and temporary files. .gridss.working subdirectories will be created. Defaults to the current directory.
	-b/--blacklist: BED file containing regions to ignore
	-s/--steps: processing steps to run. Defaults to all steps. Multiple steps are specified using comma separators. Possible steps are: setupreference, preprocess, assemble, call, all
	-c/--configuration: configuration file use to override default GRIDSS settings.
	-l/--labels: comma separated labels to use in the output VCF for the input files. Supporting read counts for input files with the same label are aggregated (useful for multiple sequencing runs of the same sample). Labels default to input filenames, unless a single read group with a non-empty sample name exists in which case the read group sample name is used (which can be disabled by \"useReadGroupSampleNameCategoryLabel=false\" in the configuration file). If labels are specified, they must be specified for all input files.
	--externalaligner: use the system version of bwa instead of the in-process version packaged with GRIDSS
	--jvmheap: size of JVM heap for assembly and variant calling. Defaults to 27.5g to ensure GRIDSS runs on all cloud instances with approximate 32gb memory including DNANexus azure:mem2_ssd1_x8.
	--maxcoverage: maximum coverage. Regions with coverage in excess of this are ignored.
	--picardoptions: additional standard Picard command line options. Useful options include VALIDATION_STRINGENCY=LENIENT and COMPRESSION_LEVEL=0. See https://broadinstitute.github.io/picard/command-line-overview.html
	--useproperpair: use SAM 'proper pair' flag to determine whether a read pair is discordant. Default: use library fragment size distribution to determine read pair concordance
	--concordantreadpairdistribution: portion of 6 sigma read pairs distribution considered concordantly mapped. Default: 0.995
	--keepTempFiles: keep intermediate files. Not recommended except for debugging due to the high disk usage.
	"

OPTIONS=r:o:a:t:j:w:b:s:c:l:
LONGOPTS=reference:,output:,assembly:,threads:,jar:,workingdir:,jvmheap:,blacklist:,steps:,configuration:,maxcoverage:,labels:,picardoptions:,jobindex:,jobnodes:,useproperpair,concordantreadpairdistribution:,keepTempFiles,sanityCheck,externalaligner
! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
	echo "$USAGE_MESSAGE" 1>&2
    exit 2
fi
eval set -- "$PARSED"
workingdir="."
reference=""
output_vcf=""
assembly=""
threads=8
jvmheap="25g"
blacklist=""
metricsrecords=10000000
steps="all"
config_file=""
maxcoverage=50000
labels=""
picardoptions=""
jobindex="0"
jobnodes="1"
useproperpair="false"
readpairpdistribution="0.995"
keepTempFiles="false"
sanityCheck="false"
externalaligner="false"
while true; do
    case "$1" in
        -r|--reference)
            reference="$2"
            shift 2
            ;;
		-l|--labels)
            labels="$2"
            shift 2
            ;;
		-s|--steps)
            steps="${2,,}"
            shift 2
            ;;
		-w|--workingdir)
            workingdir="$2"
            shift 2
            ;;
        -o|--output)
            output_vcf="$2"
            shift 2
            ;;
		-a|--assembly)
            assembly="$2"
            shift 2
            ;;
		-b|--blacklist)
            blacklist="$2"
            shift 2
            ;;
		-j|--jar)
            GRIDSS_JAR="$2"
            shift 2
            ;;
		--jvmheap)
            jvmheap="$2"
            shift 2
            ;;
		-t|--threads)
			printf -v threads '%d\n' "$2" 2>/dev/null
			printf -v threads '%d' "$2" 2>/dev/null
			shift 2
			;;
		-c|--configuration)
			config_file=$2
			shift 2
			;;
		--maxcoverage)
			maxcoverage=$2
			shift 2
			;;
		--picardoptions)
			picardoptions=$2
			shift 2
			;;
		--jobindex)
			jobindex=$2
			shift 2
			;;
		--jobnodes)
			jobnodes=$2
			shift 2
			;;
		--concordantreadpairdistribution)
			readpairpdistribution=$2
			shift 2
			;;
		--useproperpair)
			useproperpair="true"
			shift 1
			;;
		--keepTempFiles)
			keepTempFiles="true"
			shift 1
			;;
		--sanityCheck)
			sanityCheck="true"
			shift 1
			;;
		--externalaligner)
			externalaligner="true"
			shift 1
			;;
		--)
            shift
            break
            ;;
        *)
            echo "Programming error"
            exit 3
            ;;
    esac
done
do_setupreference=false
do_preprocess=false
do_assemble=false
do_call=false
for step in $(echo $steps | tr ',' ' ' ) ; do
	if [[ "$step" == "all" ]] ; then
		do_setupreference=true
		do_preprocess=true
		do_assemble=true
		do_call=true
	elif [[ "$step" == "setupreference" ]] ; then
		do_setupreference=true
	elif [[ "$step" == "preprocess" ]] ; then
		do_preprocess=true
	elif [[ "$step" == "assemble" ]] ; then
		do_assemble=true
	elif [[ "$step" == "call" ]] ; then
		do_call=true
	else
		echo "Unknown step \"$step\"" 1>&2
		exit 18
	fi
done
### Find the jars
find_jar() {
	env_name=$1
	if [[ -f "${!env_name:-}" ]] ; then
		echo "${!env_name}"
	else
		echo "Unable to find $2 jar. Specify using the environment variant $env_name, or the --jar command line parameter." 1>&2
		exit 1
	fi
}
gridss_jar=$(find_jar GRIDSS_JAR gridss)
##### --workingdir
echo "Using working directory \"$workingdir\"" 1>&2
if [[ "$workingdir" == "" ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Working directory must be specified. Specify using the --workingdir command line argument" 1>&2
	exit 3
fi
if [[ "$(tr -d ' 	\n' <<< "$workingdir")" != "$workingdir" ]] ; then
		echo "workingdir cannot contain whitespace" 1>&2
		exit 16
	fi
if [[ ! -d $workingdir ]] ; then
	mkdir -p $workingdir
	if [[ ! -d $workingdir ]] ; then
		echo Unable to create working directory $workingdir 1>&2
		exit 2
	fi
fi
workingdir=$(dirname $workingdir/placeholder)
##### --reference
echo "Using reference genome \"$reference\"" 1>&2
if [[ "$reference" == "" ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Reference genome must be specified. Specify using the --reference command line argument" 1>&2
	exit 6
fi
if [ ! -f $reference ] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Missing reference genome $reference. Specify reference location using the --reference command line argument" 1>&2
	exit 6
fi

##### --assembly
if [[ $do_assemble == "true" ]] ; then
	if [[ "$assembly" == "" ]] ; then
		echo "$USAGE_MESSAGE"  1>&2
		echo "Specify assembly bam location using the --assembly command line argument. Assembly location must be in a writeable directory." 1>&2
		exit 4
	fi
	mkdir -p $(dirname $assembly)
	if [[ ! -d $(dirname $assembly) ]] ; then
		echo Unable to parent create directory for $assembly 1>&2
		exit 2
	fi
	echo "Using assembly bam $assembly" 1>&2
fi

##### --output
if [[ $do_call == "true" ]] ; then
	if [[ "$output_vcf" == "" ]] ; then
		echo "$USAGE_MESSAGE"  1>&2
		echo "Output VCF not specified. Use --output to specify output file." 1>&2
		exit 9
	fi
	mkdir -p $(dirname $output_vcf)
	if [[ ! -d $(dirname $output_vcf) ]] ; then
		echo "Unable to create directory for $output_vcf for output VCF." 1>&2
		exit 2
	fi
	echo "Using output VCF $output_vcf" 1>&2
fi
##### --threads
if [[ "$threads" -lt 1 ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Illegal thread count: $threads. Specify an integer thread count using the --threads command line argument" 1>&2
	exit 10
fi
if [[ "$threads" -gt 8 ]] ; then
	echo "WARNING: GRIDSS scales sub-linearly at high thread count. Up to 8 threads is the recommended level of parallelism." 1>&2
fi
echo "Using $threads worker threads." 1>&2
if [[ "$blacklist" == "" ]] ; then
	blacklist_arg=""
	echo "Using no blacklist bed. The encode DAC blacklist is recommended for hg19." 1>&2
elif [[ ! -f $blacklist ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Missing blacklist file $blacklist" 1>&2
	exit 11
else
	blacklist_arg="BLACKLIST=$blacklist"
	echo "Using blacklist $blacklist" 1>&2
	if [[ "$(tr -d ' 	\n' <<< "$blacklist_arg")" != "$blacklist_arg" ]] ; then
		echo "blacklist cannot contain whitespace" 1>&2
		exit 16
	fi
fi
if [[ "$jvmheap" == "" ]] ; then
	if [[ $threads -gt 8 ]] ; then
		echo "Warning: GRIDSS assembly may stall and run out of memory. with $threads and $jvmheap heap size." 1>&2
	fi
fi
echo "Using JVM maximum heap size of $jvmheap for assembly and variant calling." 1>&2
if [[ "$@" == "" ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "At least one input bam must be specified." 1>&2
fi
for f in $@ ; do
	if [[ ! -f $f ]] ; then
		echo "Input file $f does not exist"  1>&2
		exit 12
	fi
done
config_args=""
if [[ "$config_file" != "" ]] ; then
	if [[ ! -f $config_file ]] ; then
	echo "Configuration file $config_file does not exist"  1>&2
		exit 13
	fi
	config_args="CONFIGURATION_FILE=$config_file"
fi
input_args=""
for f in $@ ; do
	if [[ "$(tr -d ' 	\n' <<< "$f")" != "$f" ]] ; then
		echo "input filenames and paths cannot contain whitespace" 1>&2
		exit 16
	fi
	echo "Using input file $f" 1>&2
	input_args="$input_args INPUT=$f"
done
if [[ "$labels" != "" ]] ; then
	nows_labels=$(tr -d ' 	\n' <<< "$labels")
	if [[ "$nows_labels" != "$labels" ]] ; then
		echo "input labels cannot contain whitespace" 1>&2
		exit 16
	fi
	IFS=',' read -ra LABEL_ARRAY  <<< "$nows_labels"
	for label in "${LABEL_ARRAY[@]}" ; do
		input_args="$input_args INPUT_LABEL=$label"
		echo label is $label
	done
	
fi

# Validate tools exist on path
for tool in Rscript samtools java ; do
	if ! which $tool >/dev/null; then echo "Error: unable to find $tool on \$PATH" 1>&2 ; exit 2; fi
	echo "Found $(which $tool)" 1>&2 
done
if which /usr/bin/time >/dev/null ; then
	timecmd="/usr/bin/time"
	echo "Found /usr/bin/time" 1>&2
else
	timecmd=""
	echo "Not found /usr/bin/time" 1>&2
fi
samtools --version-only 1>&2 || ( echo "Your samtools version is so old it does not support --version-only. Update samtools." 1>&2 && exit 19 )
Rscript --version 1>&2
if [[ "$externalaligner" == "true" ]] ; then
	if ! which bwa >/dev/null; then echo "Error: unable to find bwa on \$PATH" 1>&2 ; exit 2; fi
	echo "bwa $(bwa 2>&1 | grep Version || echo -n)" 1>&2
fi
if which /usr/bin/time >/dev/null ; then
	/usr/bin/time --version 1>&2
fi
/bin/bash --version 2>&1 | head -1 1>&2

# check java version is ok using the gridss.Echo entry point
if java -cp $gridss_jar gridss.Echo ; then
	java -version 1>&2
else
	echo "Unable to run GRIDSS jar. GRIDSS requires java 1.8 or later." 1>&2
	java -version 1>&2
	exit 14
fi

if ! java -Xms$jvmheap -cp $gridss_jar gridss.Echo ; then
	echo "Failure invoking java with --jvmheap parameter of \"$jvmheap\". Specify a JVM heap size (e.g. \"31g\") that is valid for this machine."
	exit 15
fi

timestamp=$(date +%Y%m%d_%H%M%S)
logfile=$workingdir/gridss.full.$timestamp.$HOSTNAME.$$.log
timinglogfile=$workingdir/gridss.timing.$timestamp.$HOSTNAME.$$.log
if [[ "$timecmd" != "" ]] ; then
	timecmd="/usr/bin/time --verbose -a -o $timinglogfile"
	if ! $timecmd echo 2>&1 > /dev/null; then
		timecmd="/usr/bin/time -a -o $timinglogfile"
	fi
	if ! $timecmd echo 2>&1 > /dev/null ; then
		timecmd=""
		echo "Unexpected /usr/bin/time version. Not logging timing information." 1>&2
	fi
fi

ulimit -n $(ulimit -Hn) # Reduce likelihood of running out of open file handles 
unset DISPLAY # Prevents errors attempting to connecting to an X server when starting the R plotting device
echo "Max file handles: $(ulimit -n)" 1>&2 

echo -n "$(date)	Running GRIDSS steps" 1>&2
if [[ $do_setupreference == "true" ]] ; then
	echo -n " setupreference," 1>&2
fi
if [[ $do_preprocess == "true" ]] ; then
	echo -n " preprocess," 1>&2
fi
if [[ $do_assemble == "true" ]] ; then
	echo -n " assemble," 1>&2
fi
if [[ $do_call == "true" ]] ; then
	echo -n " call," 1>&2
fi
echo ". The full log is in $logfile" 2>&1 

# For debugging purposes, we want to keep all our 
if [[ $keepTempFiles == "true" ]] ; then
	rmcmd="echo rm"
	jvm_args="-Dgridss.keepTempFiles=true"
else
	rmcmd="rm"
	jvm_args=""
fi

jvm_args="$jvm_args \
	-Dreference_fasta=$reference \
	-Dsamjdk.use_async_io_read_samtools=true \
	-Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dsamjdk.buffer_size=4194304"

readpairing_args=""

READ_PAIR_CONCORDANT_PERCENT="$readpairpdistribution"
if [[ "$useproperpair" == "true" ]] ; then
	readpairing_args="READ_PAIR_CONCORDANT_PERCENT=null"
fi

if [[ $do_setupreference == true ]] ; then
	if [[ ! -f ${reference}.fai ]] && [[ ! -f $(basename $reference .fa).fai ]] && [[ ! -f $(basename $reference .fasta).fai ]]  ; then
		echo "$(date)	samtools faidx	(once-off setup for reference genome)" | tee -a $timinglogfile
		$timecmd samtools faidx $reference 1>&2 2>> $logfile
	fi
	if [[ ! -f ${reference}.bwa  ]] ; then
		echo "$(date)	bwa index	(once-off setup for reference genome)" | tee -a $timinglogfile
		$timecmd bwa index $reference 1>&2 2>> $logfile
	fi
	echo "$(date)	PrepareReference	(once-off setup for reference genome)" | tee -a $timinglogfile
	$timecmd java -Xmx4g $jvm_args \
		-cp $gridss_jar gridss.PrepareReference \
		REFERENCE_SEQUENCE=$reference \
		1>&2 2>> $logfile
fi
exit

if [[ $do_preprocess == true ]] ; then
	if [[ "$jobnodes" != "1" ]] ; then
		echo "Error: Preprocessing does not support multiple nodes for a given input file." 2>&1
		echo "	To perform parallel per-input preprocessing, run independent preprocessing per input file." 2>&1
		exit 16
	fi
	for f in $@ ; do
		echo "$(date)	Start pre-processing	$f"
		dir=$workingdir/$(basename $f).gridss.working
		prefix=$workingdir/$(basename $f).gridss.working/$(basename $f)
		tmp_prefix=$workingdir/$(basename $f).gridss.working/tmp.$(basename $f)
		mkdir -p $dir
		if [[ ! -d $dir ]] ; then
			echo Unable to create directory $dir 1>&2
			exit 2
		fi
		if [[ ! -f $prefix.sv.bam ]] ; then
			echo "$(date)	CollectInsertSizeMetrics	$f	first $metricsrecords records" | tee -a $timinglogfile
			{ $timecmd java -Xmx4g $jvm_args \
					-cp $gridss_jar gridss.analysis.CollectGridssMetrics \
					TMP_DIR=$dir \
					ASSUME_SORTED=true \
					I=$f \
					O=$tmp_prefix \
					THRESHOLD_COVERAGE=$maxcoverage \
					FILE_EXTENSION=null \
					GRIDSS_PROGRAM=null \
					PROGRAM=null \
					PROGRAM=CollectInsertSizeMetrics \
					STOP_AFTER=$metricsrecords \
					$picardoptions \
			; } 1>&2 2>> $logfile
			echo "$(date)	CollectGridssMetricsAndExtractSVReads|samtools	$f" | tee -a $timinglogfile
			{ $timecmd java -Xmx4g $jvm_args \
					-cp $gridss_jar gridss.CollectGridssMetricsAndExtractSVReads \
					TMP_DIR=$dir \
					ASSUME_SORTED=true \
					I=$f \
					O=$prefix \
					THRESHOLD_COVERAGE=$maxcoverage \
					FILE_EXTENSION=null \
					GRIDSS_PROGRAM=null \
					GRIDSS_PROGRAM=CollectCigarMetrics \
					GRIDSS_PROGRAM=CollectMapqMetrics \
					GRIDSS_PROGRAM=CollectTagMetrics \
					GRIDSS_PROGRAM=CollectIdsvMetrics \
					GRIDSS_PROGRAM=ReportThresholdCoverage \
					PROGRAM=null \
					PROGRAM=CollectInsertSizeMetrics \
					SV_OUTPUT=/dev/stdout \
					COMPRESSION_LEVEL=0 \
					METRICS_OUTPUT=$prefix.sv_metrics \
					INSERT_SIZE_METRICS=$tmp_prefix.insert_size_metrics \
					$readpairing_args \
					UNMAPPED_READS=false \
					MIN_CLIP_LENGTH=5 \
					INCLUDE_DUPLICATES=true \
					$picardoptions \
			| $timecmd samtools sort \
					-n \
					-T $tmp_prefix.namedsorted-tmp \
					-Obam \
					-o $tmp_prefix.namedsorted.bam \
					-@ $threads \
					/dev/stdin \
			; } 1>&2 2>> $logfile
			if [[ -f $tmp_prefix.insert_size_metrics ]] ; then
				$rmcmd $tmp_prefix.insert_size_metrics $tmp_prefix.insert_size_histogram.pdf
			fi
			if [[ "$externalaligner" == "true" ]] ; then
				echo "$(date)	ComputeSamTags|samtools	$f" | tee -a $timinglogfile
				{ $timecmd java -Xmx4g $jvm_args \
						-cp $gridss_jar gridss.ComputeSamTags \
						TMP_DIR=$dir \
						WORKING_DIR=$workingdir \
						REFERENCE_SEQUENCE=$reference \
						COMPRESSION_LEVEL=0 \
						I=$tmp_prefix.namedsorted.bam \
						O=/dev/stdout \
						RECALCULATE_SA_SUPPLEMENTARY=true \
						SOFTEN_HARD_CLIPS=true \
						FIX_MATE_INFORMATION=true \
						FIX_DUPLICATE_FLAG=true \
						FIX_SA=true \
						FIX_MISSING_HARD_CLIP=true \
						TAGS=null \
						TAGS=NM \
						TAGS=SA \
						TAGS=R2 \
						TAGS=Q2 \
						TAGS=MC \
						TAGS=MQ \
						ASSUME_SORTED=true \
						$picardoptions \
				| $timecmd samtools sort \
						-T $tmp_prefix.coordinate-tmp \
						-Obam \
						-o $tmp_prefix.coordinate.bam \
						-@ $threads \
						/dev/stdin \
				; } 1>&2 2>> $logfile
				$rmcmd $tmp_prefix.namedsorted.bam
				echo "$(date)	SoftClipsToSplitReads	$f" | tee -a $timinglogfile
				{ $timecmd java -Xmx4g $jvm_args \
						-Dsamjdk.create_index=false \
						-Dgridss.gridss.output_to_temp_file=true \
						-cp $gridss_jar gridss.SoftClipsToSplitReads \
						TMP_DIR=$workingdir \
						WORKING_DIR=$workingdir \
						REFERENCE_SEQUENCE=$reference \
						I=$tmp_prefix.coordinate.bam \
						O=$tmp_prefix.sc2sr.primary.sv.bam \
						OUTPUT_UNORDERED_RECORDS=$tmp_prefix.sc2sr.supp.sv.bam \
						WORKER_THREADS=$threads \
						$picardoptions \
				&& $rmcmd $tmp_prefix.coordinate.bam \
				&& $timecmd samtools sort \
						-@ $threads \
						-T $tmp_prefix.sc2sr.suppsorted.sv-tmp \
						-Obam \
						-o $tmp_prefix.sc2sr.suppsorted.sv.bam \
						$tmp_prefix.sc2sr.supp.sv.bam \
				&& $rmcmd $tmp_prefix.sc2sr.supp.sv.bam \
				&& $timecmd samtools merge \
						-@ $threads \
						$prefix.sv.tmp.bam \
						$tmp_prefix.sc2sr.primary.sv.bam \
						$tmp_prefix.sc2sr.suppsorted.sv.bam \
				&& $timecmd samtools index $prefix.sv.tmp.bam \
				&& $rmcmd $tmp_prefix.sc2sr.primary.sv.bam \
				&& $rmcmd $tmp_prefix.sc2sr.suppsorted.sv.bam \
				&& mv $prefix.sv.tmp.bam $prefix.sv.bam \
				&& mv $prefix.sv.tmp.bam.bai $prefix.sv.bam.bai \
				; } 1>&2 2>> $logfile
			else
				echo "$(date)	PreprocessForBreakendAssembly|samtools	$f" | tee -a $timinglogfile
				{ $timecmd java -Xmx4g $jvm_args \
						-cp $gridss_jar gridss.PreprocessForBreakendAssembly \
						TMP_DIR=$dir \
						WORKING_DIR=$workingdir \
						REFERENCE_SEQUENCE=$reference \
						COMPRESSION_LEVEL=0 \
						I=$tmp_prefix.namedsorted.bam \
						O=/dev/stdout \
						WORKER_THREADS=$threads \
						ALIGNER=BWAMEM \
						ALIGNER_BATCH_SIZE=10000 \
						$picardoptions \
				&& $timecmd samtools sort \
						-@ $threads \
						-T $tmp_prefix.sc2sr.suppsorted.sv-tmp \
						-Obam \
						-o $prefix.sv.tmp.bam \
						/dev/stdin \
				&& $timecmd samtools index $prefix.sv.tmp.bam \
				&& mv $prefix.sv.tmp.bam $prefix.sv.bam \
				&& mv $prefix.sv.tmp.bam.bai $prefix.sv.bam.bai \
				; } 1>&2 2>> $logfile
			fi
			echo "$(date)	Complete pre-processing	$f"
		else
			echo "$(date)	Skipping pre-processing as $prefix.sv.bam already exists. $f"
		fi
	done
else
	echo "$(date)	Skipping pre-processing."
fi
if [[ $sanityCheck == "true" ]] ; then 
	echo "$(date)	Sanity checking *.sv.bam"
	java -Xmx$jvmheap $jvm_args \
		-cp $gridss_jar gridss.SanityCheckEvidence \
		TMP_DIR=$workingdir \
		WORKING_DIR=$workingdir \
		REFERENCE_SEQUENCE=$reference \
		WORKER_THREADS=$threads \
		$input_args \
		$blacklist_arg \
		$config_args \
		ASSEMBLY=ignored \
		OUTPUT_ERROR_READ_NAMES=reads_failing_sanity_check.txt \
		1>&2 2>> $logfile
fi
if [[ $do_assemble == true ]] ; then
	echo "$(date)	Start assembly	$assembly" | tee -a $timinglogfile
	if [[ ! -f $assembly ]] ; then
		echo "$(date)	AssembleBreakends	$assembly	job $jobindex	total jobs $jobnodes" | tee -a $timinglogfile
		{ $timecmd java -Xmx$jvmheap $jvm_args \
				-Dgridss.gridss.output_to_temp_file=true \
				-cp $gridss_jar gridss.AssembleBreakends \
				JOB_INDEX=$jobindex \
				JOB_NODES=$jobnodes \
				TMP_DIR=$workingdir \
				WORKING_DIR=$workingdir \
				REFERENCE_SEQUENCE=$reference \
				WORKER_THREADS=$threads \
				O=$assembly \
				$input_args \
				$blacklist_arg \
				$config_args \
				$picardoptions \
				$readpairing_args \
		; } 1>&2 2>> $logfile
	else
		echo "$(date)	Skipping assembly as $assembly already exists.	$assembly"
	fi
	if [[ "$jobnodes" != "1" ]] ; then
		echo "Assembly processing for job index $jobindex complete." 2>&1 
		echo "To perform the gather/reduce after all jobs are complete, run assembly without specifying --jobnodes." 2>&1 
		exit 0
	fi
	dir=$workingdir/$(basename $assembly).gridss.working/
	prefix=$dir/$(basename $assembly)
	tmp_prefix=$dir/tmp.$(basename $assembly)
	if [[ ! -f $prefix.sv.bam ]] ; then
		echo "$(date)	CollectGridssMetrics	$assembly" | tee -a $timinglogfile
		{ $timecmd java -Xmx4g $jvm_args \
				-cp $gridss_jar gridss.analysis.CollectGridssMetrics \
				I=$assembly \
				O=$prefix \
				REFERENCE_SEQUENCE=$reference \
				THRESHOLD_COVERAGE=$maxcoverage \
				TMP_DIR=$workingdir \
				FILE_EXTENSION=null \
				GRIDSS_PROGRAM=null \
				GRIDSS_PROGRAM=CollectCigarMetrics \
				GRIDSS_PROGRAM=CollectMapqMetrics \
				GRIDSS_PROGRAM=CollectTagMetrics \
				GRIDSS_PROGRAM=CollectIdsvMetrics \
				GRIDSS_PROGRAM=ReportThresholdCoverage \
				PROGRAM=null \
				PROGRAM=QualityScoreDistribution \
				$picardoptions \
		; } 1>&2 2>> $logfile
		if [[ "$externalaligner" == "true" ]] ; then
			echo "$(date)	SoftClipsToSplitReads	$assembly" | tee -a $timinglogfile
			{ $timecmd java -Xmx4g $jvm_args \
					-Dgridss.async.buffersize=16 \
					-Dsamjdk.create_index=false \
					-Dgridss.gridss.output_to_temp_file=true \
					-cp $gridss_jar gridss.SoftClipsToSplitReads \
					TMP_DIR=$dir \
					WORKING_DIR=$workingdir \
					REFERENCE_SEQUENCE=$reference \
					WORKER_THREADS=$threads \
					I=$assembly \
					O=$tmp_prefix.sc2sr.primary.sv.bam \
					OUTPUT_UNORDERED_RECORDS=$tmp_prefix.sc2sr.supp.sv.bam \
					READJUST_PRIMARY_ALIGNMENT_POSITON=true \
					$picardoptions \
			&& $timecmd samtools sort \
					-@ $threads \
					-T $tmp_prefix.sc2sr.suppsorted.sv-tmp \
					-Obam \
					-o $tmp_prefix.sc2sr.suppsorted.sv.bam \
					$tmp_prefix.sc2sr.supp.sv.bam \
			&& $rmcmd $tmp_prefix.sc2sr.supp.sv.bam \
			&& $timecmd samtools merge \
					-@ $threads \
					$prefix.sv.tmp.bam \
					$tmp_prefix.sc2sr.primary.sv.bam \
					$tmp_prefix.sc2sr.suppsorted.sv.bam \
			&& $timecmd samtools index $prefix.sv.tmp.bam \
			&& $rmcmd $tmp_prefix.sc2sr.primary.sv.bam \
			&& $rmcmd $tmp_prefix.sc2sr.suppsorted.sv.bam \
			&& mv $prefix.sv.tmp.bam $prefix.sv.bam \
			&& mv $prefix.sv.tmp.bam.bai $prefix.sv.bam.bai \
			; } 1>&2 2>> $logfile
		else
			{ $timecmd java -Xmx4g $jvm_args \
					-Dgridss.async.buffersize=16 \
					-Dsamjdk.create_index=false \
					-Dgridss.gridss.output_to_temp_file=true \
					-cp $gridss_jar gridss.SoftClipsToSplitReads \
					TMP_DIR=$dir \
					WORKING_DIR=$workingdir \
					REFERENCE_SEQUENCE=$reference \
					WORKER_THREADS=$threads \
					I=$assembly \
					O=/dev/stdout \
					ALIGNER=BWAMEM \
					ALIGNER_BATCH_SIZE=10000 \
					REALIGN_ANCHORING_BASES=true \
					READJUST_PRIMARY_ALIGNMENT_POSITON=true \
					COMPRESSION_LEVEL=0 \
					$picardoptions \
			&& $timecmd samtools sort \
					-@ $threads \
					-T $tmp_prefix.sc2sr.suppsorted.sv-tmp \
					-Obam \
					-o $prefix.sv.tmp.bam \
					/dev/stdin \
			&& $timecmd samtools index $prefix.sv.tmp.bam \
			&& mv $prefix.sv.tmp.bam $prefix.sv.bam \
			&& mv $prefix.sv.tmp.bam.bai $prefix.sv.bam.bai \
			; } 1>&2 2>> $logfile
		fi
	fi
	echo "$(date)	Complete assembly	$assembly"
else
	echo "$(date)	Skipping assembly	$assembly"
fi
if [[ $sanityCheck == "true" ]] ; then 
	java -Xmx$jvmheap $jvm_args \
		-cp $gridss_jar gridss.SanityCheckEvidence \
		TMP_DIR=$workingdir \
		WORKING_DIR=$workingdir \
		REFERENCE_SEQUENCE=$reference \
		WORKER_THREADS=$threads \
		$input_args \
		$blacklist_arg \
		$config_args \
		ASSEMBLY=$assembly \
		$readpairing_args \
		OUTPUT_ERROR_READ_NAMES=reads_failing_sanity_check.txt
fi
if [[ $do_call == true ]] ; then
	echo "$(date)	Start calling	$output_vcf" | tee -a $timinglogfile
	if [[ "$jobnodes" != "1" ]] ; then
		echo "Error: variant calling does not (yet) support multiple nodes for a given input file." 2>&1
		exit 16
	fi
	if [[ ! -f $output_vcf ]] ; then
		dir=$workingdir/$(basename $output_vcf).gridss.working
		prefix=$dir/$(basename $output_vcf)
		mkdir -p $dir
		if [[ ! -d $dir ]] ; then
			echo Unable to create directory $dir 1>&2
			exit 2
		fi
		echo "$(date)	IdentifyVariants	$output_vcf" | tee -a $timinglogfile
		{ $timecmd java -Xmx$jvmheap $jvm_args \
				-Dgridss.output_to_temp_file=true \
				-cp $gridss_jar gridss.IdentifyVariants \
				TMP_DIR=$workingdir \
				WORKING_DIR=$workingdir \
				REFERENCE_SEQUENCE=$reference \
				WORKER_THREADS=$threads \
				$input_args \
				$blacklist_arg \
				$config_args \
				ASSEMBLY=$assembly \
				OUTPUT_VCF=$prefix.unallocated.vcf \
				$readpairing_args \
		; } 1>&2 2>> $logfile
		echo "$(date)	AnnotateVariants	$output_vcf" | tee -a $timinglogfile
		{ $timecmd java -Xmx$jvmheap $jvm_args \
				-Dgridss.output_to_temp_file=true \
				-cp $gridss_jar gridss.AnnotateVariants \
				TMP_DIR=$workingdir \
				WORKING_DIR=$workingdir \
				REFERENCE_SEQUENCE=$reference \
				WORKER_THREADS=$threads \
				$input_args \
				$blacklist_arg \
				$config_args \
				ASSEMBLY=$assembly \
				INPUT_VCF=$prefix.unallocated.vcf \
				OUTPUT_VCF=$prefix.allocated.vcf \
				$picardoptions \
				$readpairing_args \
		; } 1>&2 2>> $logfile
		$rmcmd $prefix.unallocated.vcf
		echo "$(date)	AnnotateUntemplatedSequence	$output_vcf" | tee -a $timinglogfile
		{ $timecmd java -Xmx4g $jvm_args \
				-Dgridss.output_to_temp_file=true \
				-cp $gridss_jar gridss.AnnotateUntemplatedSequence \
				TMP_DIR=$workingdir \
				WORKING_DIR=$workingdir \
				REFERENCE_SEQUENCE=$reference \
				WORKER_THREADS=$threads \
				INPUT=$prefix.allocated.vcf \
				OUTPUT=$output_vcf \
				$picardoptions \
		; } 1>&2 2>> $logfile
		$rmcmd $prefix.allocated.vcf
	else
		echo "$(date)	Skipping variant calling	$output_vcf"
	fi
	echo "$(date)	Complete calling	$output_vcf" | tee -a $timinglogfile
fi
if [[ -f $logfile ]] ; then
	echo "$(date)	Run complete with $(grep WARNING $logfile | wc -l) warnings and $(grep ERROR $logfile | wc -l) errors."
else 
	echo "$(date)	Run complete with no log file. Output file already exists?"
fi
trap - EXIT
exit 0 # success!