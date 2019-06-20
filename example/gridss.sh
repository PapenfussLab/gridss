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
USAGE_MESSAGE="gridss.sh --reference <reference.fa> --output <output.vcf> --assembly <assembly.bam> [--threads n] [--jar gridss.jar] [--workingdir <directory>] [--jvmheap <threads * 4>g] [--blacklist <blacklist.bed>] input1.bam [input2.bam [...]]"

OPTIONS=r:o:a:t:j:w:b:
LONGOPTS=reference:,output:,assembly:,threads:,jar:,workingdir:,jvmheap:,blacklist:
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
threads=$(nproc)
gridss_jar=""
jvmheap=""
blacklist=""
while true; do
    case "$1" in
        -r|--reference)
            reference="$2"
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
if [[ ! -f $gridss_jar ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Unable to find GRIDSS jar. Specify location using the --jar command line argument" 1>&2
	exit 2
fi
gridss_jar=$(readlink -f $gridss_jar || echo -n)
echo "Using GRIDSS jar $gridss_jar" 1>&2
if [[ "$workingdir" == "" ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Working directory must be specified. Specify using the --workingdir command line argument" 1>&2
	exit 2
else
	if [[ ! -d $workingdir ]] ; then
		if ! mkdir -p $workingdir ; then
			echo Unable to create working directory $workingdir 1>&2
			exit 2
		fi
	fi
	workingdir=$(readlink -f $workingdir || echo -n)
	echo "Using working directory $workingdir" 1>&2
fi
if [[ "$assembly" == "" ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Specify assembly bam location using the --assembly command line argument. Assembly location must be in a writeable directory." 1>&2
fi
assembly=$(readlink -f $assembly || echo -n)
if [[ "$reference" == "" ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Specify reference location using the --reference command line argument" 1>&2
fi
reference=$(readlink -f $reference || echo -n)
if [[ ! -f "$reference" ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Missing reference genome $reference" 1>&2
fi
if [[ ! -f ${reference}.fai ]] && [[ ! -f ${reference/.fa/.fai} ]] && [[ ! -f ${reference/.faasta/.fai} ]]  ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Unable to find fai index for reference genome." 1>&2
	echo "Please create using `samtools faidx $reference`" 1>&2
fi
if [[ ! -f ${reference}.bwt  ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Unable to find bwa index ${reference}.bwt for reference genome." 1>&2
	echo "Please create using `bwa index $reference`" 1>&2
fi
echo "Using reference genome $reference" 1>&2
output_vcf=$(readlink -f $output_vcf || echo -n)
if [[ "$output_vcf" == "" ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Output VCF must be specified in a writable directory. Specify using the --output command line argument" 1>&2
	exit 2
fi 
echo "Using output VCF $output_vcf" 1>&2
if [[ "$threads" -lt 1 ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Illegal thread count: $threads. Specify an integer thread count using the --threads command line argument" 1>&2
	exit 2
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
	exit 2
else
	blacklist_arg="BLACKLIST=$(readlink -f $blacklist)"
	echo "Using blacklist $blacklist" 1>&2
fi
if [[ "$jvmheap" == "" ]] ; then
	if [[ $threads -gt 7 ]] ; then
		jvmheap="31g"
	else 
		jvmheap="$((threads * 4))g"
	fi
fi
echo "Using JVM maximum heap size of $jvmheap" 1>&2
if [[ "$@" == "" ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "At least one input bam must be specified." 1>&2
fi
input_files=""
for F in $@ ; do
	F=$(readlink -f $F)
	input_files="$input_files $F"
done
input_args=""
for F in $input_files ; do
	echo "Using input file $F" 1>&2
	input_args="$input_args INPUT=$F"
done

# Validate tools exist on path
if ! which bwa >/dev/null; then echo "Error: unable to find bwa on \$PATH" ; exit 2; fi
if ! which java >/dev/null; then echo "Error: unable to find java on \$PATH" ; exit 2; fi
if ! which Rscript >/dev/null; then echo "Error: unable to find Rscript on \$PATH" ; exit 2; fi
if ! which /usr/bin/time >/dev/null; then echo "Error: unable to find /usr/bin/time" ; exit 2; fi

mkdir -p $workingdir
logfile=$workingdir/gridss.lite.$HOSTNAME.$$.log
timinglogfile=$workingdir/gridss.timing.$HOSTNAME.$$.log

ulimit -n $(ulimit -Hn) # Reduce likelihood of running out of open file handles 
unset DISPLAY # Prevents errors attempting to connecting to an X server when starting the R plotting device

# -Dreference_fasta is only required for CRAM input files
# -Dgridss.gridss.output_to_temp_file=true allows GRIDSS to continue where it left off without data errors due to truncated files
# -Dsamjdk.create_index=true is required for multi-threaded operation
# -Dsamjdk.use_async_io allow for async read/write on background threads which improves BAM I/O performancce
echo "gridss.analysis.CallVariants" >> $timinglogfile
/usr/bin/time -a -o $timinglogfile java -Xmx$jvmheap \
	-Dreference_fasta="$reference" \
	-Dsamjdk.create_index=true \
	-Dsamjdk.use_async_io_read_samtools=true \
	-Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dsamjdk.buffer_size=4194304 \
	-Dgridss.gridss.output_to_temp_file=true \
	-cp $gridss_jar gridss.CallVariants \
	TMP_DIR=$workingdir \
	WORKING_DIR=$workingdir \
	REFERENCE_SEQUENCE="$reference" \
	$input_args \
	OUTPUT="$output_vcf" \
	ASSEMBLY="$assembly" \
	BLACKLIST="$blacklist" \
	WORKER_THREADS=$threads \
	2>&1 | tee -a $logfile

