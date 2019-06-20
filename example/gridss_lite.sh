#!/bin/bash
#
# GRIDSS lite: a high speed GRIDSS pipeline
#
# This pipeline trades slightly lower sensitivity for improved runtime performance
# by only performing assembly around putative high confidence breaks
#
set -o errexit -o pipefail -o noclobber -o nounset
! getopt --test > /dev/null 
if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
    echo '`getopt --test` failed in this environment.'
    exit 1
fi
USAGE_MESSAGE="gridss_lite.sh --reference <reference.fa> --output <output.vcf> --assembly <assembly.bam> [--threads n] [--jar gridss.jar] [--workingdir <directory>] [--jvmheap <threads * 4>g] [--blacklist <blacklist.bed>] [--firstpassqual <250>] input1.bam [input2.bam [...]]"

OPTIONS=r:o:a:t:j:w:b:
LONGOPTS=reference:,output:,assembly:,threads:,jar:,workingdir:,padding:,metricsrecords:,jvmheap:,blacklist:,firstpassqual:
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
padding=2000
metricsrecords=10000000
firstpassqual=250
cleanup="y"
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
		--padding)
			printf -v padding '%d\n' "$2" 2>/dev/null
			printf -v padding '%d' "$2" 2>/dev/null
            shift 2
            ;;
		--metricsrecords)
			printf -v metricsrecords '%d\n' "$2" 2>/dev/null
			printf -v metricsrecords '%d' "$2" 2>/dev/null
            shift 2
            ;;
		--firstpassqual)
			printf -v firstpassqual '%d\n' "$2" 2>/dev/null
			printf -v firstpassqual '%d' "$2" 2>/dev/null
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
if [[ "$gridss_jar" == "" ]] ; then
	gridss_jar=$(ls -1 ./gridss*gridss-jar-with-dependencies.jar | tail -1)
fi
if [[ ! -f $gridss_jar ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Unable to find GRIDSS jar. Specify location using the --jar command line argument" 1>&2
	exit 2
fi
gridss_jar=$(readlink -f $gridss_jar)
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
	workingdir=$(readlink -f $workingdir)
	echo "Using working directory $workingdir" 1>&2
fi
if [[ "$assembly" == "" ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Specify assembly bam location using the --assembly command line argument" 1>&2
fi
assembly=$(readlink -f $assembly)
echo Using assembly output $assembly
if [[ "$reference" == "" ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Specify reference location using the --reference command line argument" 1>&2
fi
reference=$(readlink -f $reference)
if [[ ! -f "$reference" ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Missing reference genome $reference" 1>&2
fi
if [[ ! -f ${reference}.fai ]] && [[ ! -f ${reference/.fa/.fai} ]] && [[ ! -f ${reference/.faasta/.fai} ]]  ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Unable to find fai index for reference genome." 1>&2
	echo "Please create using `samtools faidx $reference`" 1>&2
fi
echo "Using reference genome $reference" 1>&2
output_vcf=$(readlink -f $output_vcf)
if [[ "$output_vcf" == "" ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Output VCF must be specified. Specify using the --output command line argument" 1>&2
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
if [[ "$padding" -lt 1 ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Invalid region padding size $padding." 1>&2
	exit 2
fi
if [[ "$metricsrecords" -lt 1 ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Invalid number of metrics records to approximate library distributions from." 1>&2
	exit 2
fi
if [[ "$firstpassqual" -lt 1 ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Invalid first pass quality score threshold." 1>&2
	exit 2
fi
echo "Using first pass quality score threshold of $firstpassqual" 1>&2
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
if ! which samtools >/dev/null; then echo "Error: unable to find samtools on \$PATH" ; exit 2; fi
if ! which bwa >/dev/null; then echo "Error: unable to find bwa on \$PATH" ; exit 2; fi
if ! which java >/dev/null; then echo "Error: unable to find java on \$PATH" ; exit 2; fi
if ! which /usr/bin/time >/dev/null; then echo "Error: unable to find /usr/bin/time" ; exit 2; fi

workingdir=$workingdir/$(basename ${output_vcf}.gridss.lite.working)
mkdir -p $workingdir
cd $workingdir

logfile=$workingdir/gridss.lite.$HOSTNAME.$$.log
timinglogfile=$workingdir/gridss.timing.$HOSTNAME.$$.log
firstpassvcf=$workingdir/$(basename $output_vcf).firstpass.vcf
regionbed=$workingdir/$(basename $output_vcf).firstpass.bed
jvmargs="-ea \
	-Dreference_fasta="$reference" \
	-Dsamjdk.create_index=true \
	-Dsamjdk.use_async_io_read_samtools=true \
	-Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dsamjdk.buffer_size=4194304 \
	-Dgridss.gridss.output_to_temp_file=true \
	-cp $gridss_jar "
for input in $input_files ; do
	input_working_dir=$workingdir/$(basename $input).gridss.working
	echo "Calculating metrics for $input based on first $metricsrecords reads" 1>&2
	mkdir -p $input_working_dir
	echo "gridss.analysis.CollectGridssMetrics (first pass) $input" >> $timinglogfile 
	/usr/bin/time -a -o $timinglogfile java -Xmx2g $jvmargs gridss.analysis.CollectGridssMetrics \
			TMP_DIR=$input_working_dir \
			ASSUME_SORTED=true \
			I=$input \
			O=$input_working_dir/$(basename $input) \
			THRESHOLD_COVERAGE=50000 \
			FILE_EXTENSION=null \
			GRIDSS_PROGRAM=null \
			GRIDSS_PROGRAM=CollectCigarMetrics \
			GRIDSS_PROGRAM=CollectMapqMetrics \
			GRIDSS_PROGRAM=CollectTagMetrics \
			GRIDSS_PROGRAM=CollectIdsvMetrics \
			GRIDSS_PROGRAM=ReportThresholdCoverage \
			PROGRAM=null \
			PROGRAM=CollectInsertSizeMetrics \
			STOP_AFTER=$metricsrecords \
			2>&1 | tee -a $logfile
	echo "gridss.analysis.CollectStructuralVariantReadMetrics $input" >> $timinglogfile 
	/usr/bin/time -a -o $timinglogfile java -Xmx2g $jvmargs gridss.analysis.CollectStructuralVariantReadMetrics \
			TMP_DIR=$input_working_dir \
			I=$input \
			OUTPUT=$input_working_dir/$(basename $input).sv_metrics \
			INSERT_SIZE_METRICS=$input_working_dir/$(basename $input).insert_size_metrics \
			STOP_AFTER=$metricsrecords \
			2>&1 | tee -a $logfile
	if ! grep -e "^SA" $input_working_dir/$(basename $input).tag_metrics 2> /dev/null ; then
		echo "Input file is missing SA tag. "  1>&2
		echo "GRIDSS lite requires the input file to contain split read alignments" 1>&2
		echo "such as those output by bwa mem."  1>&2
		echo "gridss.SoftClipsToSplitReads can be used to convert soft clipped reads to split reads."  1>&2
		exit 3
	fi
	if [[ ! -f $workingdir/empty.bam ]] ; then
		samtools view -H $input | samtools view -b - > $workingdir/empty.bam
	fi
	sv_bam=$input_working_dir/$(basename $input).sv.bam
	ln -s $input $sv_bam
	if [[ -f $input.bai ]] ; then
		ln -s $input.bai $sv_bam.bai
	elif [[ -f ${input/.bam/.bai} ]] ; then
		ln -s ${input/.bam/.bai} $sv_bam.bai
	else
		echo "Cannot file .bai index for $input" 1>&2
		exit 4
	fi
done

mkdir -p $workingdir/empty.bam.gridss.working
if [[ ! -f $workingdir/empty.bam.gridss.working/empty.bam.sv.bam ]] ; then
	ln -s $workingdir/empty.bam $workingdir/empty.bam.gridss.working/empty.bam.sv.bam
fi
samtools index $workingdir/empty.bam
samtools index $workingdir/empty.bam.gridss.working/empty.bam.sv.bam

#TODO: R2 error should only be output when performing assembly
#TODO: Hard clip error should only be output when performing assembly
#TODO: IdentifyVariants should write tmp file
#TODO: Use a configuration file to filter to Q250 instead of an awk script
#TODO: Add timing to maximal clique chunk timing (cut/paste from assembly timing)
# First pass: call only from RP and SR
echo "variantcalling.minScore=$firstpassqual" > $workingdir/firstpass.properties
echo "gridss.analysis.IdentifyVariants (first pass)" >> $timinglogfile 
/usr/bin/time -a -o $timinglogfile java -Xmx$jvmheap $jvmargs gridss.IdentifyVariants \
			TMP_DIR=$workingdir \
			WORKING_DIR=$workingdir \
			REFERENCE_SEQUENCE=$reference \
			$input_args \
			OUTPUT_VCF=$firstpassvcf \
			ASSEMBLY=empty.bam \
			WORKER_THREADS=$threads \
			CONFIGURATION_FILE=$workingdir/firstpass.properties \
			$blacklist_arg \
			2>&1 | tee -a $logfile
# Convert to bed format
grep -v "^#" < $firstpassvcf | awk '{print $1 "\t" $2 - 1 "\t" $2}' > $regionbed
input_args=""
for input in $input_files ; do
	# extract the reads flanking the interesting calls
	echo "gridss.analysis.ExtractFullReads (second pass)" >> $timinglogfile 
	/usr/bin/time -a -o $timinglogfile java -Xmx$jvmheap $jvmargs gridss.ExtractFullReads \
		B=$regionbed \
		REGION_PADDING_SIZE=$padding \
		COMPRESSION_LEVEL=0 \
		I=$input \
		O=$workingdir/$(basename $input) \
		2>&1 | tee -a $logfile
	# remove the symlink so GRIDSS will create a targeted .sv.bam
	# with appropriate tags (e.g. R2)
	rm $workingdir/$(basename $input).gridss.working/*.sv.ba*
	input_args="$input_args INPUT=$workingdir/$(basename $input)"
done
echo "gridss.analysis.CallVariants" >> $timinglogfile
/usr/bin/time -a -o $timinglogfile java -Xmx$jvmheap $jvmargs gridss.CallVariants \
	TMP_DIR=$workingdir \
	WORKING_DIR=$workingdir \
	REFERENCE_SEQUENCE=$reference \
	COMPRESSION_LEVEL=0 \
	$input_args \
	OUTPUT=$output_vcf \
	ASSEMBLY=$assembly \
	WORKER_THREADS=$threads \
	$blacklist_arg \
	2>&1 | tee -a $logfile