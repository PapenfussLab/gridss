#!/bin/bash
#
# Performs targeted extract of fragments overlapping a region of interest
#
# ../scripts/gridss_extract_overlapping_fragments.sh -t $(nproc) --targetvcf COLO829v003T.purple.sv.vcf -o out.example.targeted.bam ~/colo829/COLO829R_dedup.realigned.bam
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

SAMTOOLS_MAJOR_VERSION_REQUIREMENT="1"
SAMTOOLS_MINOR_VERSION_REQUIREMENT="10"

BEDTOOLS_MAJOR_VERSION_REQUIREMENT="2"
BEDTOOLS_MINOR_VERSION_REQUIREMENT="27"

workingdir="."
output_bam=""
threads=1
targetbed=""
targetvcf=""
targetmargin=5000
metricsrecords=10000000
metricsmaxcoverage=100000
USAGE_MESSAGE="
Usage: gridss_extract_overlapping_fragments.sh [options] --targetbed <target.bed> --targetvcf <target.vcf> [--targetmargin $targetmargin] input.bam

	Extract all alignments of reads overlapping any of the target regions.
	
	This utility is very similar to \"samtools view -L\" but ensures that if
	any alignment record overlaps a region of interest, all alignment records
	with that read name are included in the output regardless of their alignment
	location. This ensure that all records in a read pair or a chimeric alignment
	are jointly included or excluded from the output. This requires a full
	traversal of the output.
	
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
	-w/--workingdir: directory to place intermediate, temporary files and GRIDSS metrics.
			GRIDSS should be run using the same working directory or these metrics will be ignored (Default: $workingdir)
	-j/--jar: location of GRIDSS jar
	--metricsrecords: number of reads to process to generate metrics. (Default: $metricsrecords)
			Note that these are not a random subsample but the first N reads in the input file.
			This means that the sampling is biased towards telomeric sequence which,
			if insufficent number of reads are samples, will reduce the QUAL scores due to the
			higher soft clipping and discordant mapping rate in telomeric sequences.
"

OPTIONS=j:o:t:w:
LONGOPTS=jar:,output:,threads:,workingdir:,targetbed:,targetvcf:,targetmargin:,metricsrecords:
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
		-j|--jar)
			GRIDSS_JAR="$2"  # FIXME appears unused
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
# $1 is message to write
write_status() {
	echo "$(date): $1" 1>&2
}
if [[ ! -d $workingdir ]] ; then
	mkdir -p $workingdir
	if [[ ! -d $workingdir ]] ; then
		write_status "Unable to create working directory $workingdir"
		exit $EX_CANTCREAT
	fi
fi
if [[ $# -eq 0 ]]; then
	write_status "$USAGE_MESSAGE"
	write_status "Missing input file."
	exit $EX_USAGE
fi
if [[ $# -ne 1 ]]; then
	write_status "$USAGE_MESSAGE"
	write_status "Error: found multiple input files."
	exit $EX_USAGE
fi
if [[ ! -f "$1" ]] ; then
	write_status "$USAGE_MESSAGE"
	write_status "Missing input file."
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
	write_status "Using region of interest BED file $targetvcf"
else
	if [[ ! -f "$targetvcf" ]] ; then
		write_status "$USAGE_MESSAGE"
		write_status "Missing --targetvcf file $targetvcf"
		exit $EX_NOINPUT
	fi
	write_status "Using SVs in $targetvcf as regions of interest"
fi
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
write_status "Using GRIDSS jar $gridss_jar"
##### --output
if [[ "$output_bam" == "" ]] ; then
	output_bam="$1.targeted.bam"
fi
mkdir -p "$(dirname $output_bam)"
if [[ ! -d $(dirname $output_bam) ]] ; then
	write_status "Unable to create directory for $output_bam for output file."
	exit $EX_CANTCREAT
fi
write_status "Using output file $output_bam"
##### --threads
if [[ "$threads" -lt 1 ]] ; then
	write_status "$USAGE_MESSAGE"
	write_status "Illegal thread count: $threads. Specify an integer thread count using the --threads command line argument"
	exit $EX_USAGE
fi
write_status "Using $threads worker threads."
for tool in gridsstools samtools ; do
	if ! which $tool >/dev/null; then
		write_status "Error: unable to find $tool on \$PATH"
		exit $EX_CONFIG
	fi
	write_status "Found $(which $tool)"
done

##### samtools version check
samtools_version=$(samtools --version | grep samtools | cut -b 10-)
write_status "samtools version: $samtools_version"
samtools_major_version=$(echo "$samtools_version" | cut -f 1 -d ".")
samtools_minor_version=$(echo "$samtools_version" | cut -f 2 -d ".")
if [[ "$samtools_major_version" == "$SAMTOOLS_MAJOR_VERSION_REQUIREMENT" ]] ; then
	if [[ "$samtools_minor_version" -lt "$SAMTOOLS_MINOR_VERSION_REQUIREMENT" ]] ; then
		write_status "samtools $SAMTOOLS_MAJOR_VERSION_REQUIREMENT.$SAMTOOLS_MINOR_VERSION_REQUIREMENT or later is required"
		exit $EX_CONFIG
	fi
elif [[ "$samtools_major_version" -gt "$SAMTOOLS_MAJOR_VERSION_REQUIREMENT" ]] ; then
	# newer samtools - assume it's fine
	echo -n 
else
	write_status "samtools $SAMTOOLS_MAJOR_VERSION_REQUIREMENT.$SAMTOOLS_MINOR_VERSION_REQUIREMENT or later is required"
	exit $EX_CONFIG
fi

##### bedtools version check
# bedtools --version yields 'bedtools vX.Y.Z'
bedtools_version="$(bedtools --version | cut -d' ' -f2)"
write_status "bedtools version: $bedtools_version"
bedtools_major_version="$(echo "$bedtools_version" | cut -f1 -d"." | sed 's/v//')"
bedtools_minor_version="$(echo "$bedtools_version" | cut -f2 -d".")"

if [[ "$bedtools_major_version" == "$BEDTOOLS_MAJOR_VERSION_REQUIREMENT" ]] ; then
	if [[ "$bedtools_minor_version" -lt "$BEDTOOLS_MINOR_VERSION_REQUIREMENT" ]] ; then
		write_status "bedtools $BEDTOOLS_MAJOR_VERSION_REQUIREMENT.$BEDTOOLS_MINOR_VERSION_REQUIREMENT or later is required"
		exit $EX_CONFIG
	fi
elif [[ "$bedtools_major_version" -gt "$BEDTOOLS_MAJOR_VERSION_REQUIREMENT" ]] ; then
	# newer bedtools - assume it's fine
	echo -n
else
	write_status "bedtools $BEDTOOLS_MAJOR_VERSION_REQUIREMENT.$BEDTOOLS_MINOR_VERSION_REQUIREMENT or later is required"
	exit $EX_CONFIG
fi

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
working_prefix=$workingdir/tmp.$(basename "$output_bam").gridss
target_no_slop_file=$working_prefix.target.no_slop.bed
target_file=$working_prefix.target.bed
genome_file=$working_prefix.bedtools.genome
script_vcf_to_bed=$working_prefix.vcf2bed.R
script_bam_header_to_bed_genome="$working_prefix.bam_header_to_genome.R"
rm -f $working_prefix*
if [[ "$targetbed" == "" ]] ; then
	if which Rscript > /dev/null ; then
		write_status "Converting SV VCF to BED of breakpoint positions."
	else 
		write_status "Unable to convert SV VCF to BED. Rscript is not on PATH. "
		write_status "Note that the StructuralVariantAnnotation BioConductor package is also required."
		exit $EX_USAGE
	fi
	rm -f $script_vcf_to_bed
	cat > $script_vcf_to_bed << EOF
suppressPackageStartupMessages(library(StructuralVariantAnnotation))
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
	Rscript "$script_vcf_to_bed"
else
	cp $targetbed $target_no_slop_file
fi

# Create bedtools genome file from samtools header
# Just the SN and LN attributes in tab-delimited format
# @SQ\tSN:<seq_name>\tLN:<seq_length> ==> <seq_name>\t<seq_length>
write_status "Generating bedtools genome file from bam header"
cat > "$script_bam_header_to_bed_genome" << EOF
#!/usr/bin/env Rscript

## Get contigs from bam header into bedtools genome file format

# Imports
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(tidyverse))

# Create df / render / print
Rsamtools::scanBamHeader("$input_bam", what="targets")[["$input_bam"]][["targets"]] %>%
  # Convert to standard data frame
  as.data.frame() %>%
  # Convert rowname to column
  tibble::rownames_to_column(var="chr") %>%
  # Convert to tibble
  tibble::as_tibble() %>%
  # Format
  readr::format_delim(delim="\t", col_names=FALSE) %>%
  # And print
  cat()
EOF

# Start by pulling out the bam header
Rscript "$script_bam_header_to_bed_genome" > "$genome_file"

# Use bedtools slop to extend the bed file by margins set by targetmargin variable
bedtools slop \
  -i "$target_no_slop_file" \
  -g "$genome_file" \
  -b "$targetmargin" > "$target_file"

# Then use gridsstools to extract the reads of interest
write_status "Extracting reads of interest by $targetmargin bp"
gridsstools extractFragmentsToBam -@ $threads -o $output_bam <(samtools view -M -@ $threads -L $target_file $input_bam | cut -f 1) $input_bam
gridss_dir=$workingdir/$(basename $output_bam).gridss.working
gridss_prefix=$gridss_dir/$(basename $output_bam)
if [[ ! -f $gridss_prefix.insert_size_metrics ]] ; then
	write_status "Generating GRIDSS metrics from first $metricsrecords records"
	mkdir -p "$gridss_dir"
	jvm_args=" \
		-Dpicard.useLegacyParser=false \
		-Dsamjdk.use_async_io_read_samtools=true \
		-Dsamjdk.use_async_io_write_samtools=true \
		-Dsamjdk.use_async_io_write_tribble=true \
		-Dsamjdk.buffer_size=4194304 \
		-Dsamjdk.async_io_read_threads=$threads"
	java -Xmx2g $jvm_args \
		-cp $gridss_jar gridss.analysis.CollectGridssMetrics \
		--INPUT $input_bam \
		--OUTPUT $gridss_prefix \
		--THRESHOLD_COVERAGE $metricsmaxcoverage \
		--TMP_DIR $workingdir \
		--FILE_EXTENSION null \
		--STOP_AFTER $metricsrecords
else 
	write_status "Skipped GRIDSS metrics generation. Found	$gridss_prefix.insert_size_metrics"
fi
write_status "Done"
rm -f $working_prefix*
trap - EXIT
exit 0 # success!
