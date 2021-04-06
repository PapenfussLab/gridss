#!/bin/bash
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

output="/dev/stdout"
threads=$(nproc)
workingdir="."
rm="RepeatMasker"
rmargs="--species human"
minlength="20"
#fields="INSRM INSRMRT INSRMRC INSRMRO INSRMP"
USAGE_MESSAGE="
Usage: gridss_annotate_vcf_repeatmasker.sh [options] input.vcf
	-o,--output: output vcf file. Defaults to stdout.
	-j/--jar: location of GRIDSS jar
	-t/--threads: number of threads to use. Defaults to the number of cores available ($threads)
	-w/--workingdir: directory to place intermediate and temporary files (Default: $workingdir)
	--rm: RepeatMasker executable. (Default: $rm)
	--rmargs: additional RepeatMasker arguments (Default: $rmargs)
	--minlength: minimum length of inserted sequence to annotate. (Default: $minlength)
	"
#--fields: INFO fields to populate. (Default: $fields)
OPTIONS=o:t:j:w:
LONGOPTS=output:,workingdir:,threads:,jar:,rm:,rmargs:,minlength:
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
		-o|--output)
			output="$2"
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
		-w|--workingdir)
			workingdir="$2"
			shift 2
			;;
		--minlength)
			printf -v minlength '%d\n' "$2" 2>/dev/null
			printf -v minlength '%d' "$2" 2>/dev/null
			shift 2
			;;
		--rm)
			rm=$2
			shift 2
			;;
		--rmargs)
			rmargs=$2
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
write_status() {
	echo "$(date): $1" 1>&2
}
if [[ "$#" != "1" ]] ; then
	echo "$USAGE_MESSAGE"
	exit $EX_USAGE
fi
if [[ "$1" == "" ]] ; then
	echo "$USAGE_MESSAGE"
	write_status "Missing input vcf" 
	exit $EX_USAGE
fi
if [[ ! -f $1 ]] ; then
	
	echo "$USAGE_MESSAGE" 1>&2
	write_status "Input file '$1' not found." 1>&2
	exit $EX_NOINPUT
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
##### --threads
if [[ "$threads" -lt 1 ]] ; then
	write_status "$USAGE_MESSAGE"
	write_status "Illegal thread count: $threads. Specify an integer thread count using the --threads command line argument"
	exit $EX_USAGE
fi
write_status  "Using $threads worker threads."
if [[ "$output" == "" ]] ; then
	echo "$USAGE_MESSAGE"
	write_status "Missing output vcf. Specify with --output"
	exit $EX_USAGE
fi
# Validate tools exist on path
for tool in $rm java ; do
	if ! which $tool >/dev/null; then
		write_status "Error: unable to find $tool on \$PATH"
		exit $EX_CONFIG
	fi
	write_status "Found $(which $tool)"
done
write_status "RepeatMasker version: $($rm -v)"
write_status "bash version: $(/bin/bash --version 2>&1 | head -1)"

# check java version is ok using the gridss.Echo entry point
if java -cp $gridss_jar gridss.Echo ; then
	write_status "java version: $(java -version 2>&1)"
else
	write_status "Unable to run GRIDSS jar. GRIDSS requires java 1.8 or later."
	write_status "java version: $(java -version  2>&1)"
	exit $EX_CONFIG
fi

input_vcf="$1"
if [[ ! -f "$input_vcf" ]] ; then
	echo "$USAGE_MESSAGE"
	write_status "Missing input vcf \"$input_vcf\""
	exit $EX_NOINPUT
fi

jvm_args=" \
	-Dpicard.useLegacyParser=false \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dsamjdk.buffer_size=4194304 \
	-Dsamjdk.async_io_read_threads=$threads"

temp_fa="$workingdir/$(basename $input_vcf).fa"
mkdir -p "$workingdir"
# Extract inserted sequences
java -Xmx64m $jvm_args -cp $GRIDSS_JAR \
	gridss.InsertedSequencesToFasta \
	--INPUT "$input_vcf" \
	--OUTPUT "$temp_fa" \
	--MIN_SEQUENCE_LENGTH $minlength
if [[ ! -s $temp_fa ]] ; then
	# No inserted sequences in the VCF
	cat $input_vcf > $output
else
	# Run RepeatMasker on sequences
	$rm $rmargs -pa $threads -dir "$workingdir" "$temp_fa"
	# We can annotate with full alignment information if it exists
	rmfile="$temp_fa.align"
	if [[ ! -f $temp_fa.align ]] ; then
		rmfile="$temp_fa.out"
	fi
	# Parse and annotate
	# Needs the extra memory since the RM output needs to be loaded into memory
	java -Xmx2g $jvm_args -cp $GRIDSS_JAR \
		gridss.repeatmasker.AnnotateVariantsRepeatMasker \
		--INPUT "$input_vcf" \
		--OUTPUT "$output" \
		--REPEAT_MASKER "$rmfile"
	# --TAGS $fields
fi

trap - EXIT
exit 0 # success!