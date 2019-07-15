#!/bin/bash
#
# GRIDSS: a sensitive structural variant calling toolkit
#
# This script is a simple wrapper around the all-in-one gridss.CallVariants entry point
#
# Example ./gridss.sh  -t 1 -b wgEncodeDacMapabilityConsensusExcludable.bed -r ../hg19.fa -o gridss.full.chr12.1527326.DEL1024.vcf -a gridss.full.chr12.1527326.DEL1024.assembly.bam -j ../target/gridss-2.4.0-gridss-jar-with-dependencies.jar chr12.1527326.DEL1024.bam
set -o errexit -o pipefail -o noclobber -o nounset
! getopt --test > /dev/null 
if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
    echo '`getopt --test` failed in this environment.'
    exit 1
fi
USAGE_MESSAGE="
Usage: gridss.sh --reference <reference.fa> --output <output.vcf> --assembly <assembly.bam> [--threads n] [--jar gridss.jar] [--workingdir <directory>] [--jvmheap 28g] [--blacklist <blacklist.bed>] [--steps All|PreProcess|Assemble|Call] [--configuration gridss.properties] [--maxcoverage 50000] input1.bam [input2.bam [...]]

	-r/--reference: reference genome to use. Must have a .fai index file and a bwa index.
	-o/--output: output VCF.
	-a/--assembly: location of the GRIDSS assembly BAM. This file will be created by GRIDSS.
	-t/--threads: number of threads to use. Defaults to the number of cores available.
	-j/--jar: location of GRIDSS jar
	-w/--workingdir: directory to place GRIDSS intermediate and temporary files. .gridss.working subdirectories will be created. Defaults to the current directory.
	-b/--blacklist: BED file containing regions to ignore
	-s/--steps: processing steps to run. Defaults to all steps.
	-c/--configuration: configuration file use to override default GRIDSS settings.
	--jvmheap: size of JVM heap for assembly and variant calling. Defaults to 27.5g to ensure GRIDSS runs on all cloud instances with approximate 32gb memory including DNANexus azure:mem2_ssd1_x8.
	--maxcoverage: maximum coverage. Regions with coverage in excess of this are ignored.
	"
	

OPTIONS=r:o:a:t:j:w:b:s:c:
LONGOPTS=reference:,output:,assembly:,threads:,jar:,workingdir:,jvmheap:,blacklist:,steps:,configuration:,maxcoverage:
! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
	echo "$USAGE_MESSAGE" 1>&2
    exit 2
fi
eval set -- "$PARSED"
workingdir="./"
reference=""
output_vcf=""
assembly=""
threads=$(nproc)
gridss_jar=""
jvmheap="28g"
blacklist=""
metricsrecords=10000000
steps="All"
config_file=""
maxcoverage=50000
while true; do
    case "$1" in
        -r|--reference)
            reference="$2"
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
            gridss_jar="$2"
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
##### --jar
if [[ ! -f $gridss_jar ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Unable to find GRIDSS jar. Specify location using the --jar command line argument" 1>&2
	exit 2
fi
echo "Using GRIDSS jar $gridss_jar" 1>&2
##### --workingdir
if [[ "$workingdir" == "" ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Working directory must be specified. Specify using the --workingdir command line argument" 1>&2
	exit 3
fi
if [[ ! -d $workingdir ]] ; then
	if ! mkdir -p $workingdir ; then
		echo Unable to create working directory $workingdir 1>&2
		exit 2
	fi
fi
echo "Using working directory $workingdir" 1>&2
##### --reference
if [[ ! -f "$reference" ]] ; then
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
	echo "Please create using `bwa index $reference`" 1>&2
	exit 8
fi
echo "Using reference genome $reference" 1>&2

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
	echo "Using input file $f" 1>&2
	input_args="$input_args INPUT=$f"
done

# Validate tools exist on path
for tool in bwa Rscript /usr/bin/time sambamba java ; do
	if ! which $tool >/dev/null; then echo "Error: unable to find $tool on \$PATH" 1>&2 ; exit 2; fi
	echo "Found $(which $tool)" 1>&2 
done

# check java version is ok by testing for GRIDSS usage message
if java -cp $gridss_jar gridss.Echo 2>&1 | grep "USAGE:" >/dev/null ; then
	java -version 2>&1
else
	echo "Unable to run GRIDSS jar. GRIDSS requires java 1.8 or later." 2>&1
	java -version 2>&1
	exit 14
fi


timestamp=$(date +%Y%m%d_%H%M%S)
mkdir -p $workingdir
logfile=$workingdir/gridss.full.$timestamp.$HOSTNAME.$$.log
timinglogfile=$workingdir/gridss.timing.$timestamp.$HOSTNAME.$$.log

ulimit -n $(ulimit -Hn) # Reduce likelihood of running out of open file handles 
unset DISPLAY # Prevents errors attempting to connecting to an X server when starting the R plotting device
echo "Max file handles: $(ulimit -n)" 1>&2 

COMMON_JVM_ARGS="-Dreference_fasta=$reference \
	-Dsamjdk.create_index=true \
	-Dsamjdk.use_async_io_read_samtools=true \
	-Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dsamjdk.buffer_size=4194304 \
	-Dgridss.gridss.output_to_temp_file=true"
	
for f in $@ ; do
	echo "Pre-processing $F"
	dir=$workingdir/$(basename $f).gridss.working
	prefix=$workingdir/$(basename $f).gridss.working/$(basename $f)
	tmp_prefix=$workingdir/$(basename $f).gridss.working/tmp.$(basename $f)
	if ! mkdir -p $dir ; then
		echo Unable to create directory $dir 1>&2
		exit 2
	fi
	if [[ ! -f $prefix.sv.bam ]] ; then
		echo "$(date)	Running	CollectInsertSizeMetrics	$f	$metricsrecords records" >> $timinglogfile
		/usr/bin/time -a -o $timinglogfile java -Xmx256m $JVM_ARGS gridss.analysis.CollectGridssMetrics \
			ASSUME_SORTED=true \
			I=$f \
			O=$tmp_prefix \
			THRESHOLD_COVERAGE=$maxcoverage \
			FILE_EXTENSION=null \
			GRIDSS_PROGRAM=null \
			PROGRAM=null \
			PROGRAM=CollectInsertSizeMetrics \
			STOP_AFTER=$metricsrecords 2>&1 | tee -a $logfile
		exit
		echo "$(date)	Running	CollectGridssMetricsAndExtractSVReads|sambamba	$f" >> $timinglogfile
		/usr/bin/time -a -o $timinglogfile java -Xmx512m $JVM_ARGS gridss.analysis.CollectGridssMetricsAndExtractSVReads \
			TMP_DIR=$gridss_dir \
			ASSUME_SORTED=true \
			I=$bam \
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
			INSERT_SIZE_METRICS=$tmp_prefix \
			UNMAPPED_READS=false \
			MIN_CLIP_LENGTH=5 \
			INCLUDE_DUPLICATES=true | \
		sambamba sort --sort-picard --tmpdir $dir -m 8GB -o $tmp_prefix.coordinate.bam -t $threads - 2>&1 | tee -a $logfile
	fi
done

# -Dreference_fasta is only required for CRAM input files
# -Dgridss.gridss.output_to_temp_file=true allows GRIDSS to continue where it left off without data errors due to truncated files
# -Dsamjdk.create_index=true is required for multi-threaded operation
# -Dsamjdk.use_async_io allow for async read/write on background threads which improves BAM I/O performancce
echo "gridss.analysis.CallVariants" >> $timinglogfile
/usr/bin/time -a -o $timinglogfile java -Xmx$jvmheap \
	 \
	-cp $gridss_jar gridss.CallVariants \
	TMP_DIR=$workingdir \
	WORKING_DIR=$workingdir \
	REFERENCE_SEQUENCE="$reference" \
	$input_args \
	OUTPUT="$output_vcf" \
	ASSEMBLY="$assembly" \
	 $blacklist_arg  \
	WORKER_THREADS=$threads \
	2>&1 | tee -a $logfile

