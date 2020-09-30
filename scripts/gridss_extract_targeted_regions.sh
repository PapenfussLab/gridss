#!/bin/bash
#
# Performs targeted preprocessing for GRIDSS
#
# ../scripts/gridss_extract_targeted_regions.sh -t $(nproc) --targetvcf COLO829v003T.purple.sv.vcf -o out.example.targeted.bam ~/colo829/COLO829R_dedup.realigned.bam
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

workingdir="."
reference=""
output_bam=""
threads=1
targetbed=""
targetvcf=""
targetmargin=2000
metricsrecords=1000000
USAGE_MESSAGE="
Usage: gridss_create_targeted_input.sh [options] --targetbed <target.bed> --targetvcf <target.vcf> [--targetmargin $targetmargin] input.bam

	Creates a BAM file and associated GRIDSS intermediate files that can be
	used for targeted GRIDSS calling.
	
	If the input BED/VCF only includes one side of a breakpoint, the reference
	coverage annotations will be incorrect on the breakend outside the region
	of interest. Similarly, partially targeted Okazaki fragments will be have
	incomplete calls. These issues can be rectified by iteratively rerunning
	targeted GRIDSS on the GRIDSS output.

	--targetbed: BED regions of interest
	--targetvcf: SV VCF contains calls of interest
		If a SV VCF is used as input. Rscript must be on path and the
		StructuralVariantAnnotation package must be installed.
		All SV breakend are considered targets regardless of VCF notation
			- start and end of <DEL> <DUP> <INS> <INV> symbolic alleles,
			- start and end of all indels
			- all <BND> postions
			- SNVs and CNV records are ignored
		The targeted region of breakpoints with microhomology and imprecise
		extends to the bounds defined by CIPOS and CIEND.
	--targetmargin: size of flanking region to include around the targeted region (Default: $targetmargin)
	-o/--output: output BAM. Defaults to adding a .targeted.bam suffix to the input filename
	-t/--threads: number of threads to use. (Default: $threads)
	-j/--jar: location of GRIDSS jar
	-w/--workingdir: directory to place GRIDSS intermediate and temporary files.
		.gridss.working subdirectories will be created.
		This working directory must match the working directly used to run GRIDSS
		(Default: $workingdir)
	--metricsrecords: number of reads to estimate input file metrics from (Default: $metricsrecords)
"

OPTIONS=r:o:t:j:w:
LONGOPTS=reference:,output:,threads:,jar:,workingdir:,targetbed:,targetvcf:,targetmargin:,metricsrecords:
! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
	echo "$USAGE_MESSAGE" 1>&2
    exit $EX_USAGE
fi
eval set -- "$PARSED"
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
            output_bam="$2"
            shift 2
            ;;
		-j|--jar)
            GRIDSS_JAR="$2"
            shift 2
            ;;
		-t|--threads)
			printf -v threads '%d\n' "$2" 2>/dev/null
			printf -v threads '%d' "$2" 2>/dev/null
			shift 2
			;;
		--targetbed)
			targetbed="$2"
			shift 2
			;;
		--targetvcf)
			targetvcf="$2"
			shift 2
			;;
		--targetmargin)
			printf -v targetmargin '%d\n' "$2" 2>/dev/null
			printf -v targetmargin '%d' "$2" 2>/dev/null
			shift 2
			;;
		--metricsrecords)
			printf -v metricsrecords '%d\n' "$2" 2>/dev/null
			printf -v metricsrecords '%d' "$2" 2>/dev/null
			shift 2
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
logfile=$workingdir/gridss_create_targeted_input.full.$timestamp.$HOSTNAME.$$.log
# $1 is message to write
write_status() {
	echo "$(date): $1" | tee -a $logfile 1>&2
}
write_status "Full log file is: $logfile"
# Timing instrumentation
timinglogfile=$workingdir/gridss_create_targeted_input.timing.$timestamp.$HOSTNAME.$$.log
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
if [[ ! -v GRIDSS_JAR ]] || [[ ! -f  "$GRIDSS_JAR" ]] ; then
	write_status "Unable to find GRIDSS jar. Specify using the environment variant GRIDSS_JAR, or the --jar command line parameter."
	exit $EX_NOINPUT
fi
if [[ $# -eq 0 ]]; then
	write_status "$USAGE_MESSAGE"
	write_status "Missing input bam."
	exit $EX_USAGE
fi
if [[ $# -ne 1 ]]; then
	write_status "$USAGE_MESSAGE"
	write_status "Error: found multiple input files."
	exit $EX_USAGE
fi
if [[ ! -f "$1" ]] ; then
	write_status "$USAGE_MESSAGE"
	write_status "Missing input bam."
	exit $EX_NOINPUT
fi
input_bam="$1"
write_status "Using input file \"$input_bam\""
if [[ "$targetbed" != "" ]] ; then
	if [[ "$targetvcf" != "" ]] ; then
		write_status "$USAGE_MESSAGE"
		write_status "--targetbed and --targetvcf are mutually exclusive. Specify only one."
		exit $EX_USAGE
	fi
	if [[ ! -f "$targetbed" ]] ; then
		write_status "$USAGE_MESSAGE"
		write_status "Missing --targetbed file $targetbed"
		exit $EX_NOINPUT
	fi
else
	if [[ ! -f "$targetvcf" ]] ; then
		write_status "$USAGE_MESSAGE"
		write_status "Missing --targetvcf file $targetvcf"
		exit $EX_NOINPUT
	fi
fi
##### --output
if [[ "$output_bam" == "" ]] ; then
	output_bam="$1.targeted.bam"
fi
mkdir -p $(dirname $output_bam)
if [[ ! -d $(dirname $output_bam) ]] ; then
	write_status "Unable to create directory for $output_bam for output BAM."
	exit $EX_CANTCREAT
fi
write_status "Using output BAM $output_bam"
##### --threads
if [[ "$threads" -lt 1 ]] ; then
	write_status "$USAGE_MESSAGE"
	write_status "Illegal thread count: $threads. Specify an integer thread count using the --threads command line argument"
	exit $EX_USAGE
fi
write_status "Using $threads worker threads."

jvm_args="-ea \
	-Dpicard.useLegacyParser=false \
	-Dsamjdk.create_index=true \
	-Dsamjdk.use_async_io_read_samtools=true \
	-Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dsamjdk.buffer_size=4194304 \
	-Dgridss.output_to_temp_file=true \
	-cp $GRIDSS_JAR "
if [[ -f "$reference" ]] ; then
	jvm_args="$jvm_args -Dsamjdk.reference_fasta=$reference"
fi
gridss_working_dir=$workingdir/$(basename "$output_bam").gridss.working
gridss_working_prefix=$gridss_working_dir/$(basename "$output_bam")
target_no_slop_file=$gridss_working_prefix.target.no_slop.bed
target_file=$gridss_working_prefix.target.bed
script_vcf_to_bed=$gridss_working_prefix.vcf2bed.R
mkdir -p $gridss_working_dir
if [[ "$targetbed" == "" ]] ; then
	if which Rscript > /dev/null ; then
		write_status "Converting SV VCF to BED format."
	else 
		write_status "Unable to convert SV VCF to BED. Rscript is not on PATH."
		exit $EX_USAGE
	fi
	rm -f $script_vcf_to_bed
	cat > $script_vcf_to_bed << EOF
library(StructuralVariantAnnotation)
vcf = readVcf("$targetvcf")
bpgr = breakpointRanges(vcf, unpartneredBreakends=FALSE, inferMissingBreakends=TRUE)
begr = breakpointRanges(vcf, unpartneredBreakends=TRUE)
remove(vcf)
gr = c(bpgr, begr)
remove(begr)
remove(bpgr)
gr = sort(gr, ignore.strand=TRUE)
export(object=gr, con="$target_no_slop_file", format="bed")
EOF
	{ $timecmd Rscript "$script_vcf_to_bed" ; } 1>&2 2>> $logfile
else
	cp $targetbed $target_no_slop_file
fi
grep -v "^#" $target_no_slop_file | grep -v "^browser" | grep -v "^track" | awk "{OFS=\"\t\" print \$1,\$2-$targetmargin,\$3+$targetmargin}" > $target_file

if [[ ! -f $gridss_working_dir/$(basename "$output_bam").sv_metrics ]] ; then
	write_status "Calculating metrics based on first $metricsrecords reads."
	{ $timecmd java -Xmx3500m $jvm_args gridss.analysis.CollectGridssMetrics \
		--TMP_DIR "$gridss_working_dir "\
		--ASSUME_SORTED true \
		--I "$input_bam" \
		--O "$gridss_working_prefix" \
		--THRESHOLD_COVERAGE 1000000000 \
		--FILE_EXTENSION null \
		--GRIDSS_PROGRAM null \
		--GRIDSS_PROGRAM CollectCigarMetrics \
		--GRIDSS_PROGRAM CollectMapqMetrics \
		--GRIDSS_PROGRAM CollectTagMetrics \
		--GRIDSS_PROGRAM CollectIdsvMetrics \
		--GRIDSS_PROGRAM ReportThresholdCoverage \
		--PROGRAM null \
		--PROGRAM CollectInsertSizeMetrics \
		--STOP_AFTER $metricsrecords \
	&& java -Xmx8g $jvm_args gridss.analysis.CollectStructuralVariantReadMetrics \
			--TMP_DIR "$gridss_working_dir" \
			--I "$input_bam" \
			--OUTPUT "$gridss_working_prefix.sv_metrics" \
			--INSERT_SIZE_METRICS "$gridss_working_prefix.insert_size_metrics" \
			--STOP_AFTER $metricsrecords \
	; } 1>&2 2>> $logfile
else
	write_status "Skipping calculate metrics step."
fi
write_status "Extracting fragments overlapping targeted regions	$output_bam"
{ $timecmd java -Xmx3500m $jvm_args \
	gridss.ExtractFullReads \
	--B "$target_file" \
	--REGION_PADDING_SIZE $targetmargin \
	--I "$input_bam" \
	--O "$output_bam" \
; } 1>&2 2>> $logfile
write_status "To generate GRIDSS results, run GRIDSS on \"$output_bam\" with --workingdir \"$workingdir\""
trap - EXIT
exit 0 # success!
