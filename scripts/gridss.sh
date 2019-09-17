#!/bin/bash
#
# GRIDSS: a sensitive structural variant calling toolkit
#
# Example ./gridss.sh  -t 4 -b wgEncodeDacMapabilityConsensusExcludable.bed -r ../hg19.fa -w out -o out/gridss.full.chr12.1527326.DEL1024.vcf -a out/gridss.full.chr12.1527326.DEL1024.assembly.bam -j ../target/gridss-2.5.1-gridss-jar-with-dependencies.jar chr12.1527326.DEL1024.bam
set -o errexit -o pipefail -o noclobber -o nounset
! getopt --test > /dev/null 
if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
    echo '`getopt --test` failed in this environment.'
    exit 1
fi
USAGE_MESSAGE="
Usage: gridss.sh --reference <reference.fa> --output <output.vcf.gz> --assembly <assembly.bam> [--threads n] [--jar gridss.jar] [--workingdir <directory>] [--jvmheap 25g] [--blacklist <blacklist.bed>] [--steps All|PreProcess|Assemble|Call] [--configuration gridss.properties] [--maxcoverage 50000] [--labels input1,input2,...] input1.bam [input2.bam [...]]

	-r/--reference: reference genome to use. Must have a .fai index file and a bwa index.
	-o/--output: output VCF.
	-a/--assembly: location of the GRIDSS assembly BAM. This file will be created by GRIDSS.
	-t/--threads: number of threads to use. Defaults to the number of cores available.
	-j/--jar: location of GRIDSS jar
	-w/--workingdir: directory to place GRIDSS intermediate and temporary files. .gridss.working subdirectories will be created. Defaults to the current directory.
	-b/--blacklist: BED file containing regions to ignore
	-s/--steps: processing steps to run. Defaults to all steps. Multiple steps are specified using comma separators
	-c/--configuration: configuration file use to override default GRIDSS settings.
	-l/--labels: comma separated labels to use in the output VCF for the input files. Supporting read counts for input files with the same label are aggregated (useful for multiple sequencing runs of the same sample). Labels default to input filenames, unless a single read group with a non-empty sample name exists in which case the read group sample name is used (which can be disabled by \"useReadGroupSampleNameCategoryLabel=false\" in the configuration file). If labels are specified, they must be specified for all input files.
	--jvmheap: size of JVM heap for assembly and variant calling. Defaults to 27.5g to ensure GRIDSS runs on all cloud instances with approximate 32gb memory including DNANexus azure:mem2_ssd1_x8.
	--maxcoverage: maximum coverage. Regions with coverage in excess of this are ignored.
	"
	

OPTIONS=r:o:a:t:j:w:b:s:c:l:
LONGOPTS=reference:,output:,assembly:,threads:,jar:,workingdir:,jvmheap:,blacklist:,steps:,configuration:,maxcoverage:,labels:
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
threads=$(nproc)
jvmheap="25g"
blacklist=""
metricsrecords=10000000
steps="all"
config_file=""
maxcoverage=50000
labels=""
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
do_preprocess=false
do_assemble=false
do_call=false
if [[ "$steps" == *"all"* ]] ; then
	do_preprocess=true
	do_assemble=true
	do_call=true
fi
if [[ "$steps" == *"preprocess"* ]] ; then
	do_preprocess=true
fi
if [[ "$steps" == *"assemble"* ]] ; then
	do_assemble=true
fi
if [[ "$steps" == *"call"* ]] ; then
	do_call=true
fi
### Find the jars
find_jar() {
	env_name=$1
	if [[ -f "${!env_name:-}" ]] ; then
		echo "${!env_name}"
	else
		echo "Unable to find $2 jar. Specify using the environment variant $env_name" 1>&2
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
	if ! mkdir -p $workingdir ; then
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
if [[ ! -f ${reference}.fai ]] && [[ ! -f $(basename $reference .fa).fai ]] && [[ ! -f $(basename $reference .fasta).fai ]]  ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Unable to find fai index for reference genome." 1>&2
	echo "Please create using `samtools faidx $reference`" 1>&2
	exit 7
fi
if [[ ! -f ${reference}.bwt  ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Unable to find bwa index ${reference}.bwt for reference genome." 1>&2
	echo "Please create using bwa index $reference" 1>&2
	exit 8
fi

##### --assembly
if [[ $do_assemble == "true" ]] ; then
	if [[ "$assembly" == "" ]] ; then
		echo "$USAGE_MESSAGE"  1>&2
		echo "Specify assembly bam location using the --assembly command line argument. Assembly location must be in a writeable directory." 1>&2
		exit 4
	fi
	mkdir -p $(dirname $assembly) || echo "Unable to create directory $(dirname $assembly) for assembly BAM." 1>&2
	echo "Using assembly bam $assembly" 1>&2
fi

##### --output
if [[ $do_call == "true" ]] ; then
	if [[ "$output_vcf" == "" ]] ; then
		echo "$USAGE_MESSAGE"  1>&2
		echo "Output VCF not specified. Use --output to specify output file." 1>&2
		exit 9
	fi
	mkdir -p $(dirname $output_vcf) || echo "Unable to create directory $(dirname $output_vcf) for output VCF." 1>&2
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
for tool in bwa Rscript /usr/bin/time sambamba java ; do
	if ! which $tool >/dev/null; then echo "Error: unable to find $tool on \$PATH" 1>&2 ; exit 2; fi
	echo "Found $(which $tool)" 1>&2 
done
sambamba 2>&1 | grep sambamba | head -1 1>&2 || echo -n
Rscript --version 1>&2
echo "bwa $(bwa 2>&1 | grep Version || echo -n)" 1>&2
/usr/bin/time --version 1>&2
/bin/bash --version 2>&1 | head -1 1>&2

# check java version is ok by testing for GRIDSS usage message
if java -cp $gridss_jar gridss.Echo ; then
	java -version 1>&2
else
	echo "Unable to run GRIDSS jar. GRIDSS requires java 1.8 or later." 2>&1
	java -version 1>&2
	exit 14
fi

if ! java -Xms$jvmheap -cp $gridss_jar gridss.Echo ; then
	echo "Failure invoking java with --jvmheap parameter of \"$jvmheap\". Specify a JVM heap size (e.g. \"31g\") that is valid for this machine."
	exit 15
fi


timestamp=$(date +%Y%m%d_%H%M%S)
mkdir -p $workingdir
logfile=$workingdir/gridss.full.$timestamp.$HOSTNAME.$$.log
timinglogfile=$workingdir/gridss.timing.$timestamp.$HOSTNAME.$$.log

ulimit -n $(ulimit -Hn) # Reduce likelihood of running out of open file handles 
unset DISPLAY # Prevents errors attempting to connecting to an X server when starting the R plotting device
echo "Max file handles: $(ulimit -n)" 1>&2 

echo "$(date)	Running GRIDSS. The full log is in $logfile"

jvm_args="
	-Dreference_fasta=$reference \
	-Dsamjdk.use_async_io_read_samtools=true \
	-Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dsamjdk.buffer_size=4194304"

if [[ $do_preprocess == true ]] ; then
	for f in $@ ; do
		echo "$(date)	Start pre-processing	$f"
		dir=$workingdir/$(basename $f).gridss.working
		prefix=$workingdir/$(basename $f).gridss.working/$(basename $f)
		tmp_prefix=$workingdir/$(basename $f).gridss.working/tmp.$(basename $f)
		if ! mkdir -p $dir ; then
			echo Unable to create directory $dir 1>&2
			exit 2
		fi
		if [[ ! -f $prefix.sv.bam ]] ; then
			echo "$(date)	CollectInsertSizeMetrics	$f	first $metricsrecords records" | tee -a $timinglogfile
			{ /usr/bin/time -a -o $timinglogfile java -Xmx4g $jvm_args \
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
			; } 1>&2 2>> $logfile
			echo "$(date)	CollectGridssMetricsAndExtractSVReads|sambamba	$f" | tee -a $timinglogfile
			{ /usr/bin/time -a -o $timinglogfile \
				java -Xmx4g $jvm_args \
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
					UNMAPPED_READS=false \
					MIN_CLIP_LENGTH=5 \
					INCLUDE_DUPLICATES=true \
			| /usr/bin/time -a -o $timinglogfile \
				sambamba sort \
					-n \
					--tmpdir $dir \
					-m 8GB \
					-o $tmp_prefix.namedsorted.bam \
					-t $threads \
					/dev/stdin \
			; } 1>&2 2>> $logfile
			rm $tmp_prefix.insert_size_metrics $tmp_prefix.insert_size_histogram.pdf
			echo "$(date)	ComputeSamTags|sambamba	$f" | tee -a $timinglogfile
			{ /usr/bin/time -a -o $timinglogfile \
				java -Xmx4g $jvm_args \
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
					TAGS=null \
					TAGS=NM \
					TAGS=SA \
					TAGS=R2 \
					TAGS=Q2 \
					TAGS=MC \
					TAGS=MQ \
					ASSUME_SORTED=true \
			| /usr/bin/time -a -o $timinglogfile \
				sambamba sort \
					--tmpdir $dir \
					-m 8GB \
					-o $tmp_prefix.coordinate.bam \
					-t $threads \
					/dev/stdin \
			; } 1>&2 2>> $logfile
			rm $tmp_prefix.namedsorted.bam
			echo "$(date)	SoftClipsToSplitReads/bwa	$f" | tee -a $timinglogfile
			{ /usr/bin/time -a -o $timinglogfile \
				java -Xmx4g $jvm_args \
					-Dsamjdk.create_index=true \
					-Dgridss.gridss.output_to_temp_file=true \
					-cp $gridss_jar gridss.SoftClipsToSplitReads \
					TMP_DIR=$workingdir \
					WORKING_DIR=$workingdir \
					REFERENCE_SEQUENCE=$reference \
					I=$tmp_prefix.coordinate.bam \
					O=$prefix.sv.bam \
					WORKER_THREADS=$threads \
			; } 1>&2 2>> $logfile
			rm $tmp_prefix.coordinate.bam $tmp_prefix.coordinate.bam.bai
			echo "$(date)	Complete pre-processing	$f"
		else
			echo "$(date)	Skipping pre-processing as $prefix.sv.bam already exists. $f"
		fi
	done
else
	echo "$(date)	Skipping pre-processing."
fi

if [[ $do_assemble == true ]] ; then
	echo "$(date)	Start assembly	$assembly" | tee -a $timinglogfile
	if [[ ! -f $assembly ]] ; then
		echo "$(date)	AssembleBreakends	$assembly" | tee -a $timinglogfile
		{ /usr/bin/time -a -o $timinglogfile \
			java -Xmx$jvmheap $jvm_args \
				-Dgridss.gridss.output_to_temp_file=true \
				-cp $gridss_jar gridss.AssembleBreakends \
				TMP_DIR=$workingdir \
				WORKING_DIR=$workingdir \
				REFERENCE_SEQUENCE=$reference \
				WORKER_THREADS=$threads \
				O=$assembly \
				$input_args \
				$blacklist_arg \
				$config_args \
		; } 1>&2 2>> $logfile
	else
		echo "$(date)	Skipping assembly as $assembly already exists.	$assembly"
	fi
	prefix=$workingdir/$(basename $assembly).gridss.working/$(basename $assembly)
	if [[ ! -f $prefix.sv.bam ]] ; then
		echo "$(date)	CollectGridssMetrics	$assembly" | tee -a $timinglogfile
		{ /usr/bin/time -a -o $timinglogfile \
			java -Xmx4g $jvm_args \
				-cp $gridss_jar gridss.analysis.CollectGridssMetrics \
				I=$assembly \
				O=$prefix \
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
				PROGRAM=CollectAlignmentSummaryMetrics \
		; } 1>&2 2>> $logfile
		echo "$(date)	SoftClipsToSplitReads	$assembly" | tee -a $timinglogfile
		{ /usr/bin/time -a -o $timinglogfile \
			java -Xmx6g $jvm_args \
				-Dgridss.async.buffersize=16 \
				-Dsamjdk.create_index=true \
				-Dgridss.gridss.output_to_temp_file=true \
				-cp $gridss_jar gridss.SoftClipsToSplitReads \
				TMP_DIR=$workingdir \
				WORKING_DIR=$workingdir \
				REFERENCE_SEQUENCE=$reference \
				WORKER_THREADS=$threads \
				I=$assembly \
				O=$prefix.sv.bam \
				REALIGN_ENTIRE_READ=true \
		; } 1>&2 2>> $logfile
	fi
	echo "$(date)	Complete assembly	$assembly"
else
	echo "$(date)	Skipping assembly	$assembly"
fi
if [[ $do_call == true ]] ; then
	echo "$(date)	Start calling	$output_vcf" | tee -a $timinglogfile
	if [[ ! -f $output_vcf ]] ; then
		dir=$workingdir/$(basename $output_vcf).gridss.working
		prefix=$dir/$(basename $output_vcf)
		if ! mkdir -p $dir ; then
			echo Unable to create directory $dir 1>&2
			exit 2
		fi
		echo "$(date)	IdentifyVariants	$output_vcf" | tee -a $timinglogfile
		{ /usr/bin/time -a -o $timinglogfile \
			java -Xmx$jvmheap $jvm_args \
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
		; } 1>&2 2>> $logfile
		echo "$(date)	AnnotateVariants	$output_vcf" | tee -a $timinglogfile
		{ /usr/bin/time -a -o $timinglogfile \
			java -Xmx$jvmheap $jvm_args \
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
		; } 1>&2 2>> $logfile
		rm $prefix.unallocated.vcf
		echo "$(date)	AnnotateUntemplatedSequence	$output_vcf" | tee -a $timinglogfile
		{ /usr/bin/time -a -o $timinglogfile \
			java -Xmx4g $jvm_args \
				-Dgridss.output_to_temp_file=true \
				-cp $gridss_jar gridss.AnnotateUntemplatedSequence \
				TMP_DIR=$workingdir \
				WORKING_DIR=$workingdir \
				REFERENCE_SEQUENCE=$reference \
				WORKER_THREADS=$threads \
				INPUT=$prefix.allocated.vcf \
				OUTPUT=$output_vcf \
		; } 1>&2 2>> $logfile
		rm $prefix.allocated.vcf
	else
		echo "$(date)	Skipping variant calling	$output_vcf"
	fi
	echo "$(date)	Complete calling	$output_vcf" | tee -a $timinglogfile
fi
echo "$(date)	Run complete with $(grep WARNING $logfile | wc -l) warnings and $(grep ERROR $logfile | wc -l) errors."
