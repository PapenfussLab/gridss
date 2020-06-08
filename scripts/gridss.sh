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

EX_USAGE=64
EX_NOINPUT=66
EX_CANTCREAT=73
EX_CONFIG=78

USAGE_MESSAGE="
Usage: gridss.sh --reference <reference.fa> --output <output.vcf.gz> --assembly <assembly.bam> [--threads n] [--jar gridss.jar] [--workingdir <directory>] [--jvmheap 25g] [--blacklist <blacklist.bed>] [--steps All|PreProcess|Assemble|Call] [--configuration gridss.properties] [--maxcoverage 50000] [--labels input1,input2,...] input1.bam [input2.bam [...]]

	-r/--reference: reference genome to use. Must have a .fai index file and a bwa index.
	-o/--output: output VCF.
	-a/--assembly: location of the GRIDSS assembly BAM. This file will be created by GRIDSS.
	-t/--threads: number of threads to use. Defaults to the number of cores available.
	-j/--jar: location of GRIDSS jar
	-w/--workingdir: directory to place GRIDSS intermediate and temporary files. .gridss.working subdirectories will be created. Defaults to the current directory.
	-b/--blacklist: BED file containing regions to ignore
	--repeatmaskerbed: bedops rmsk2bed BED file for genome.
	-s/--steps: processing steps to run. Defaults to all steps. Multiple steps are specified using comma separators. Possible steps are: setupreference, preprocess, assemble, call, all. WARNING: multiple instances of GRIDSS generating reference files at the same time will result in file corruption. Make sure these files are generated before runninng parallel GRIDSS jobs.
	-c/--configuration: configuration file use to override default GRIDSS settings.
	-l/--labels: comma separated labels to use in the output VCF for the input files. Supporting read counts for input files with the same label are aggregated (useful for multiple sequencing runs of the same sample). Labels default to input filenames, unless a single read group with a non-empty sample name exists in which case the read group sample name is used (which can be disabled by \"useReadGroupSampleNameCategoryLabel=false\" in the configuration file). If labels are specified, they must be specified for all input files.
	--externalaligner: use the system version of bwa instead of the in-process version packaged with GRIDSS
	--jvmheap: size of JVM heap for assembly and variant calling. Defaults to 27.5g to ensure GRIDSS runs on all cloud instances with approximate 32gb memory including DNANexus azure:mem2_ssd1_x8.
	--maxcoverage: maximum coverage. Regions with coverage in excess of this are ignored.
	--picardoptions: additional standard Picard command line options. Useful options include VALIDATION_STRINGENCY=LENIENT and COMPRESSION_LEVEL=0. See https://broadinstitute.github.io/picard/command-line-overview.html
	--useproperpair: use SAM 'proper pair' flag to determine whether a read pair is discordant. Default: use library fragment size distribution to determine read pair concordance
	--concordantreadpairdistribution: portion of 6 sigma read pairs distribution considered concordantly mapped. Default: 0.995
	--keepTempFiles: keep intermediate files. Not recommended except for debugging due to the high disk usage.
	--nojni: do not use JNI native code acceleration libraries (snappy, GKL, ssw, bwa).
	--jobindex: zero-based assembly job index (only required when performing parallel assembly across multiple computers)
	--jobnodes: total number of assembly jobs (only required when performing parallel assembly across multiple computers). Note than an assembly jobs is required after all indexed jobs have been completed to gather the output files together.
	"

OPTIONS=r:o:a:t:j:w:b:s:c:l:
LONGOPTS=reference:,output:,assembly:,threads:,jar:,workingdir:,jvmheap:,blacklist:,steps:,configuration:,maxcoverage:,labels:,picardoptions:,jobindex:,jobnodes:,useproperpair,concordantreadpairdistribution:,keepTempFiles,sanityCheck,externalaligner,nojni,repeatmaskerbed:
! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
	echo "$USAGE_MESSAGE" 1>&2
    exit $EX_USAGE
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
nojni="false"
repeatmaskerbed=""
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
		--repeatmaskerbed)
			repeatmaskerbed=$2
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
		--nojni)
			nojni="true"
			shift 1
			;;
		--)
            shift
            break
            ;;
        *)
            echo "Programming error"
            exit 1
            ;;
    esac
done
do_setupreference=false
do_preprocess=false
do_assemble=false
do_call=false

##### --workingdir
echo "Using working directory \"$workingdir\"" 1>&2
if [[ "$workingdir" == "" ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Working directory must be specified. Specify using the --workingdir command line argument" 1>&2
	exit $EX_USAGE
fi
if [[ "$(tr -d ' 	\n' <<< "$workingdir")" != "$workingdir" ]] ; then
		echo "workingdir cannot contain whitespace" 1>&2
		exit $EX_USAGE
	fi
if [[ ! -d $workingdir ]] ; then
	mkdir -p $workingdir
	if [[ ! -d $workingdir ]] ; then
		echo Unable to create working directory $workingdir 1>&2
		exit $EX_CANTCREAT
	fi
fi
workingdir=$(dirname $workingdir/placeholder)
timestamp=$(date +%Y%m%d_%H%M%S)
# Logging
logfile=$workingdir/gridss.full.$timestamp.$HOSTNAME.$$.log
# $1 is message to write
write_status() {
	echo "$(date): $1" | tee -a $logfile 1>&2
}
write_status "Full log file is: $logfile"
# Timing instrumentation
timinglogfile=$workingdir/gridss.timing.$timestamp.$HOSTNAME.$$.log
if which /usr/bin/time >/dev/null ; then
	timecmd="/usr/bin/time"
	write_status "Found /usr/bin/time"
else
	timecmd=""
	write_status "Not found /usr/bin/time"
fi
if [[ "$timecmd" != "" ]] ; then
	timecmd="/usr/bin/time --verbose -a -o $timinglogfile"
	if ! $timecmd echo 2>&1 > /dev/null; then
		timecmd="/usr/bin/time -a -o $timinglogfile"
	fi
	if ! $timecmd echo 2>&1 > /dev/null ; then
		timecmd=""
		write_status "Unexpected /usr/bin/time version. Not logging timing information."
	fi
	# We don't need timing info of the echo
	rm -f $timinglogfile
fi

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
		write_status "Unknown step \"$step\""
		exit $EX_USAGE
	fi
done
### Find the jars
find_jar() {
	env_name=$1
	if [[ -f "${!env_name:-}" ]] ; then
		echo "${!env_name}"
	else
		write_status "Unable to find $2 jar. Specify using the environment variant $env_name, or the --jar command line parameter."
		exit $EX_NOINPUT
	fi
}
gridss_jar=$(find_jar GRIDSS_JAR gridss)
##### --reference
write_status "Using reference genome \"$reference\""
if [[ "$reference" == "" ]] ; then
	write_status "$USAGE_MESSAGE"
	write_status "Reference genome must be specified. Specify using the --reference command line argument"
	exit $EX_USAGE
fi
if [ ! -f $reference ] ; then
	write_status "$USAGE_MESSAGE"
	write_status "Missing reference genome $reference. Specify reference location using the --reference command line argument"
	exit $EX_USAGE
fi

##### --assembly
if [[ $do_assemble == "true" ]] ; then
	if [[ "$assembly" == "" ]] ; then
		write_status "$USAGE_MESSAGE"
		write_status "Specify assembly bam location using the --assembly command line argument. Assembly location must be in a writeable directory."
		exit $EX_USAGE
	fi
	mkdir -p $(dirname $assembly)
	if [[ ! -d $(dirname $assembly) ]] ; then
		write_status "Unable to parent create directory for $assembly"
		exit $EX_CANTCREAT
	fi
	write_status "Using assembly bam $assembly"
fi

##### --output
if [[ $do_call == "true" ]] ; then
	if [[ "$output_vcf" == "" ]] ; then
		write_status "$USAGE_MESSAGE"
		write_status "Output VCF not specified. Use --output to specify output file."
		exit $EX_USAGE
	fi
	mkdir -p $(dirname $output_vcf)
	if [[ ! -d $(dirname $output_vcf) ]] ; then
		write_status "Unable to create directory for $output_vcf for output VCF."
		exit $EX_CANTCREAT
	fi
	write_status "Using output VCF $output_vcf"
fi
##### --threads
if [[ "$threads" -lt 1 ]] ; then
	write_status "$USAGE_MESSAGE"
	write_status "Illegal thread count: $threads. Specify an integer thread count using the --threads command line argument"
	exit $EX_USAGE
fi
if [[ "$threads" -gt 8 ]] ; then
	write_status "WARNING: GRIDSS scales sub-linearly at high thread count. Up to 8 threads is the recommended level of parallelism."
fi
write_status  "Using $threads worker threads."
if [[ "$blacklist" == "" ]] ; then
	blacklist_arg=""
	write_status  "Using no blacklist bed. The encode DAC blacklist is recommended for hg19."
elif [[ ! -f $blacklist ]] ; then
	write_status  "$USAGE_MESSAGE"
	write_status  "Missing blacklist file $blacklist"
	exit $EX_NOINPUT
else
	blacklist_arg="BLACKLIST=$blacklist"
	write_status  "Using blacklist $blacklist"
	if [[ "$(tr -d ' 	\n' <<< "$blacklist_arg")" != "$blacklist_arg" ]] ; then
		write_status  "blacklist cannot contain whitespace"
		exit $EX_USAGE
	fi
fi
if [[ "$repeatmaskerbed" != "" ]] ; then
	if [[ ! -f $repeatmaskerbed ]] ; then
		write_status "RepeatMasker BED file ($repeatmaskerbed) not found."
		exit $EX_NOINPUT
	fi
fi
if [[ "$jvmheap" == "" ]] ; then
	if [[ $threads -gt 8 ]] ; then
		write_status "Warning: GRIDSS assembly may stall and run out of memory. with $threads and $jvmheap heap size."
	fi
fi
write_status  "Using JVM maximum heap size of $jvmheap for assembly and variant calling."
if [[ "$@" == "" ]] ; then
	write_status  "$USAGE_MESSAGE"
	write_status  "At least one input bam must be specified."
fi
for f in $@ ; do
	if [[ ! -f $f ]] ; then
		write_status "Input file $f does not exist"
		exit $EX_NOINPUT
	fi
done
config_args=""
if [[ "$config_file" != "" ]] ; then
	if [[ ! -f $config_file ]] ; then
		write_status "Configuration file $config_file does not exist"
		exit $EX_NOINPUT
	fi
	config_args="CONFIGURATION_FILE=$config_file"
fi
input_args=""
for f in $@ ; do
	if [[ "$(tr -d ' 	\n' <<< "$f")" != "$f" ]] ; then
		write_status "input filenames and paths cannot contain whitespace"
		exit $EX_USAGE
	fi
	write_status "Using input file $f"
	input_args="$input_args INPUT=$f"
done
if [[ "$labels" != "" ]] ; then
	nows_labels=$(tr -d ' 	\n' <<< "$labels")
	if [[ "$nows_labels" != "$labels" ]] ; then
		write_status "input labels cannot contain whitespace"
		exit $EX_USAGE
	fi
	IFS=',' read -ra LABEL_ARRAY  <<< "$nows_labels"
	for label in "${LABEL_ARRAY[@]}" ; do
		input_args="$input_args INPUT_LABEL=$label"
		write_status "label is $label"
	done
	
fi

for f1 in $@ ; do
	if [[ "$(basename $f1)" == "$(basename $assembly)" ]] ; then
		write_status "assembly and input files must have different filenames."
		exit $EX_USAGE
	fi
	for f2 in $@ ; do
		if [[ "$f1" != "$f2" ]] ; then
			if [[ "$(basename $f1)" == "$(basename $f2)" ]] ; then
				write_status "input files must have different filenames."
				exit $EX_USAGE
			fi
		fi
	done
done

# Validate tools exist on path
for tool in Rscript samtools java bwa ; do #minimap2
	if ! which $tool >/dev/null; then
		write_status "Error: unable to find $tool on \$PATH"
		exit $EX_CONFIG
	fi
	write_status "Found $(which $tool)"
done
if $(samtools --version-only 2>&1 >/dev/null) ; then
	write_status "samtools version: $(samtools --version-only 2>&1)"
else 
	write_status "Your samtools version is so old it does not support --version-only. Update samtools."
	exit $EX_CONFIG
fi
write_status "R version: $(Rscript --version 2>&1)"
write_status "bwa $(bwa 2>&1 | grep Version || echo -n)"
#write_status "minimap2 $(minimap2 --version)"
if which /usr/bin/time >/dev/null ; then
	write_status "time version: $(/usr/bin/time --version 2>&1)"
fi
write_status "bash version: $(/bin/bash --version 2>&1 | head -1)"

# check java version is ok using the gridss.Echo entry point
if java -cp $gridss_jar gridss.Echo ; then
	write_status "java version: $(java -version 2>&1)"
else
	write_status "Unable to run GRIDSS jar. GRIDSS requires java 1.8 or later."
	write_status "java version: $(java -version  2>&1)"
	exit $EX_CONFIG
fi

if ! java -Xms$jvmheap -cp $gridss_jar gridss.Echo ; then
	write_status "Failure invoking java with --jvmheap parameter of \"$jvmheap\". Specify a JVM heap size (e.g. \"31g\") that is valid for this machine."
	exit 1
fi

ulimit -n $(ulimit -Hn) # Reduce likelihood of running out of open file handles 
unset DISPLAY # Prevents errors attempting to connecting to an X server when starting the R plotting device
write_status "Max file handles: $(ulimit -n)" 1>&2 

steps_message="Running GRIDSS steps:"
if [[ $do_setupreference == "true" ]] ; then
	steps_message="$steps_message setupreference,"
fi
if [[ $do_preprocess == "true" ]] ; then
	steps_message="$steps_message preprocess,"
fi
if [[ $do_assemble == "true" ]] ; then
	steps_message="$steps_message assemble,"
fi
if [[ $do_call == "true" ]] ; then
	steps_message="$steps_message call,"
fi
write_status "$steps_message"

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

if [[ "$nojni" == "true" ]] ; then
	write_status "Disabling snappy, GKL, SSW, and in-process alignment"
	jvm_args="$jvm_args \
		-Dsnappy.disable=true
		-Dsamjdk.try_use_intel_deflater=false
		-Dsswjni.disable=true
		-Dgkljni.disable=true
		"
	externalaligner="true"
fi

aligner_args_bwa='
	ALIGNER_COMMAND_LINE=null
	ALIGNER_COMMAND_LINE=bwa
	ALIGNER_COMMAND_LINE=mem
	ALIGNER_COMMAND_LINE=-L
	ALIGNER_COMMAND_LINE=0,0
	ALIGNER_COMMAND_LINE=-t
	ALIGNER_COMMAND_LINE=%3$d
	ALIGNER_COMMAND_LINE=%2$s
	ALIGNER_COMMAND_LINE=%1$s
	'
aligner_args_bowtie2='
	ALIGNER_COMMAND_LINE=null
	ALIGNER_COMMAND_LINE=bowtie2
	ALIGNER_COMMAND_LINE=--threads
	ALIGNER_COMMAND_LINE=%3$d
	ALIGNER_COMMAND_LINE=--local
	ALIGNER_COMMAND_LINE=--mm
	ALIGNER_COMMAND_LINE=--reorder
	ALIGNER_COMMAND_LINE=-x
	ALIGNER_COMMAND_LINE=%2$s
	ALIGNER_COMMAND_LINE=-U
	ALIGNER_COMMAND_LINE=%1$s
	'
aligner_args_minimap2='
	ALIGNER_COMMAND_LINE=null
	ALIGNER_COMMAND_LINE=minimap2
	ALIGNER_COMMAND_LINE=%2$s.idx
	ALIGNER_COMMAND_LINE=-x
	ALIGNER_COMMAND_LINE=sr
	ALIGNER_COMMAND_LINE=-Y
	ALIGNER_COMMAND_LINE=-a
	ALIGNER_COMMAND_LINE=-t
	ALIGNER_COMMAND_LINE=%3$d
	ALIGNER_COMMAND_LINE=%1$s
	'

readpairing_args=""

READ_PAIR_CONCORDANT_PERCENT="$readpairpdistribution"
if [[ "$useproperpair" == "true" ]] ; then
	readpairing_args="READ_PAIR_CONCORDANT_PERCENT=null"
fi

if [[ $do_setupreference == true ]] ; then
	if [[ ! -f ${reference}.fai ]] && [[ ! -f $(basename $reference .fa).fai ]] && [[ ! -f $(basename $reference .fasta).fai ]]  ; then
		write_status "Running	samtools faidx	(once-off setup for reference genome)"
		$timecmd samtools faidx $reference 1>&2 2>> $logfile
	fi
	if [[ ! -f ${reference}.bwt  ]] ; then
		write_status "Running	bwa index	(once-off setup for reference genome)"
		$timecmd bwa index $reference 1>&2 2>> $logfile
	fi
	#if [[ ! -f ${reference}.idx  ]] ; then
	#	write_status "Running	minimap2 index	(once-off setup for reference genome)"
	#	$timecmd minimap2 -d ${reference}.idx ${reference} 1>&2 2>> $logfile
	#fi
	write_status "Running	PrepareReference	(once-off setup for reference genome)"
	$timecmd java -Xmx4g $jvm_args \
		-cp $gridss_jar gridss.PrepareReference \
		REFERENCE_SEQUENCE=$reference \
		1>&2 2>> $logfile
fi

if [[ $do_preprocess == true ]] ; then
	if [[ "$jobnodes" != "1" ]] ; then
		write_status "Error: Preprocessing does not support multiple nodes for a given input file."
		write_status "	To perform parallel per-input preprocessing, run independent preprocessing per input file."
		exit $EX_USAGE
	fi
	for f in $@ ; do
		write_status "Start pre-processing	$f"
		dir=$workingdir/$(basename $f).gridss.working
		prefix=$workingdir/$(basename $f).gridss.working/$(basename $f)
		tmp_prefix=$workingdir/$(basename $f).gridss.working/tmp.$(basename $f)
		mkdir -p $dir
		if [[ ! -d $dir ]] ; then
			write_status "Unable to create directory $dir"
			exit $EX_CANTCREAT
		fi
		if [[ ! -f $prefix.sv.bam ]] ; then
			write_status "Running	CollectInsertSizeMetrics	$f	first $metricsrecords records"
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
			write_status "Running	CollectGridssMetricsAndExtractSVReads|samtools	$f"
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
				write_status "Running	ComputeSamTags|samtools	$f"
				rm -f $tmp_prefix.coordinate-tmp*
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
				write_status "Running	SoftClipsToSplitReads	$f"
				{ $timecmd java -Xmx4g $jvm_args \
						-Dsamjdk.create_index=false \
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
				write_status "Running	PreprocessForBreakendAssembly|samtools	$f"
				rm -f $tmp_prefix.sc2sr.suppsorted.sv-tmp*
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
				| $timecmd samtools sort \
						-@ $threads \
						-T $tmp_prefix.sc2sr.suppsorted.sv-tmp \
						-Obam \
						-o $prefix.sv.tmp.bam \
						/dev/stdin \
				&& $rmcmd $tmp_prefix.namedsorted.bam \
				&& $timecmd samtools index $prefix.sv.tmp.bam \
				&& mv $prefix.sv.tmp.bam $prefix.sv.bam \
				&& mv $prefix.sv.tmp.bam.bai $prefix.sv.bam.bai \
				; } 1>&2 2>> $logfile
			fi
			if [[ ! -f $prefix.sv.bam ]] ; then
				write_status "pre-processing failed for $f"
				exit 1
			fi
			write_status "Complete pre-processing	$f"
		else
			write_status "Skipping pre-processing as $prefix.sv.bam already exists. $f"
		fi
	done
else
	write_status "Skipping pre-processing."
fi
if [[ $sanityCheck == "true" ]] ; then 
	write_status "Sanity checking *.sv.bam"
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
	write_status "Start assembly	$assembly"
	if [[ ! -f $assembly ]] ; then
		write_status "Running	AssembleBreakends	$assembly	job $jobindex	total jobs $jobnodes"
		{ $timecmd java -Xmx$jvmheap $jvm_args \
				-Dgridss.output_to_temp_file=true \
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
		write_status  "Skipping assembly as $assembly already exists.	$assembly"
	fi
	if [[ "$jobnodes" != "1" ]] ; then
		write_status "Assembly processing for job index $jobindex complete."
		write_status "To perform the gather/reduce after all jobs are complete, run assembly without specifying --jobnodes."
		exit 0
	fi
	dir=$workingdir/$(basename $assembly).gridss.working/
	prefix=$dir/$(basename $assembly)
	tmp_prefix=$dir/tmp.$(basename $assembly)
	if [[ ! -f $prefix.sv.bam ]] ; then
		write_status "Running	CollectGridssMetrics	$assembly"
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
		#externalaligner=true # TODO: get JNI minimap2 working
		if [[ "$externalaligner" == "true" ]] ; then
			write_status "Running	SoftClipsToSplitReads	$assembly"
			{ $timecmd java -Xmx4g $jvm_args \
					-Dgridss.async.buffersize=16 \
					-Dsamjdk.create_index=false \
					-Dgridss.output_to_temp_file=true \
					-cp $gridss_jar gridss.SoftClipsToSplitReads \
					TMP_DIR=$dir \
					WORKING_DIR=$workingdir \
					REFERENCE_SEQUENCE=$reference \
					WORKER_THREADS=$threads \
					I=$assembly \
					O=$tmp_prefix.sc2sr.primary.sv.bam \
					OUTPUT_UNORDERED_RECORDS=$tmp_prefix.sc2sr.supp.sv.bam \
					REALIGN_ENTIRE_READ=true \
					READJUST_PRIMARY_ALIGNMENT_POSITION=true \
					$aligner_args_bwa \
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
			write_status "Running	SoftClipsToSplitReads|samtools	$assembly"
			{ $timecmd java -Xmx4g $jvm_args \
					-Dgridss.async.buffersize=16 \
					-Dsamjdk.create_index=false \
					-cp $gridss_jar gridss.SoftClipsToSplitReads \
					TMP_DIR=$dir \
					WORKING_DIR=$workingdir \
					REFERENCE_SEQUENCE=$reference \
					WORKER_THREADS=$threads \
					I=$assembly \
					O=/dev/stdout \
					ALIGNER=BWAMEM \
					ALIGNER_BATCH_SIZE=100000 \
					REALIGN_ENTIRE_READ=true \
					READJUST_PRIMARY_ALIGNMENT_POSITION=true \
					COMPRESSION_LEVEL=0 \
					$picardoptions \
			| $timecmd samtools sort \
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
	write_status "Complete assembly	$assembly"
else
	write_status "Skipping assembly	$assembly"
fi
if [[ $sanityCheck == "true" ]] ; then 
	write_status "Running sanity checks"
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
	write_status "Start calling	$output_vcf"
	if [[ "$jobnodes" != "1" ]] ; then
		write_status "Error: variant calling does not (yet) support multiple nodes for a given input file."
		exit $EX_USAGE
	fi
	if [[ ! -f $output_vcf ]] ; then
		dir=$workingdir/$(basename $output_vcf).gridss.working
		prefix=$dir/$(basename $output_vcf)
		mkdir -p $dir
		if [[ ! -d $dir ]] ; then
			write_status "Unable to create directory $dir"
			exit $EX_CANTCREAT
		fi
		write_status "Running	IdentifyVariants	$output_vcf"
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
		write_status "Running	AnnotateVariants	$output_vcf"
		{ $timecmd java -Xmx$jvmheap $jvm_args \
				-Dgridss.output_to_temp_file=true \
				-Dgridss.async.buffersize=2048 \
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
		write_status "Running	AnnotateInsertedSequence	$output_vcf"
		repeatmaskerbed_cmdline=""
		if [[ "$repeatmaskerbed" != "" ]] ; then
			repeatmaskerbed_cmdline="REPEAT_MASKER_BED=$repeatmaskerbed"
		fi
		{ $timecmd java -Xmx4g $jvm_args \
				-Dgridss.output_to_temp_file=true \
				-cp $gridss_jar gridss.AnnotateInsertedSequence \
				TMP_DIR=$workingdir \
				WORKING_DIR=$workingdir \
				REFERENCE_SEQUENCE=$reference \
				WORKER_THREADS=$threads \
				INPUT=$prefix.allocated.vcf \
				OUTPUT=$output_vcf \
				$repeatmaskerbed_cmdline \
				$picardoptions \
		&& $rmcmd $prefix.allocated.vcf \
		; } 1>&2 2>> $logfile
	else
		write_status  "Skipping variant calling	$output_vcf"
	fi
	write_status "Complete calling	$output_vcf"
fi
if [[ -f $logfile ]] ; then
	write_status "Run complete with $(grep WARNING $logfile | wc -l) warnings and $(grep ERROR $logfile | wc -l) errors."
fi
trap - EXIT
exit 0 # success!