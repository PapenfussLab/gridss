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

kraken2db=""
workingdir="."
reference=""
output_vcf=""
threads=8
kraken2args=""
gridssargs=""
nodesdmp=""
virushostdb=""
minreads="50"
viralgenomes="1"
metricsrecords=10000000
metricsmaxcoverage=100000
maxcoverage=1000000
hosttaxid=0
force="false"
USAGE_MESSAGE="
Viral Integration Recognition Using Single Breakends

Usage: virusbreakend.sh [options] input.bam

	-r/--reference: reference genome of host species.
	-o/--output: output VCF.
	-j/--jar: location of GRIDSS jar
	-t/--threads: number of threads to use. (Default: $threads).
	-w/--workingdir: directory to place intermediate and temporary files. (Default: $workingdir).
	--kraken2db: kraken2 database
	--hosttaxid: NCBI taxonomy id of host. Used to filter viral sequences of interest to those infecting this host. Default: 0 (no host filtering)
	--virushostdb: location of virushostdb.tsv. Available from ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv
	--kraken2args: additional kraken2 arguments
	--gridssargs: additional GRIDSS arguments
	--nodesdmp: location of NCBI nodes.dmp. Can be downloaded from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip. (Default: {kraken2db}/taxonomy/nodes.dmp)
	--minreads: minimum number of viral reads perform integration detection (Default: $minreads)
	--viralgenomes: number of viral genomes to consider. Multiple closely related genomes will result in a high false negative rate due to multi-mapping reads. (Default: $viralgenomes)
	"
OPTIONS=o:t:j:w:r:f
LONGOPTS=output:,jar:,threads:,workingdir:,kraken2db:,kraken2args:,gridssargs:,nodesdmp:,minreads:,hosttaxid:,virushostdb:,force
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
			output_vcf="$2"
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
		--kraken2db)
			kraken2db="$2"
			shift 2
			;;
		--kraken2args)
			kraken2args=$2
			shift 2
			;;
		--gridssargs)
			gridssargs=$2
			shift 2
			;;
		--nodesdmp)
			nodesdmp="$2"
			shift 2
			;;
		--minreads)
			minreads="$2"
			shift 2
			;;
		--hosttaxid)
			printf -v hosttaxid '%d\n' "$2" 2>/dev/null
			printf -v hosttaxid '%d' "$2" 2>/dev/null
			shift 2
			;;
		--virushostdb)
			virushostdb="$2"
			shift 2
			;;
		-f|--force)
			force="true"
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
write_status() { # Before logging initialised
	echo "$(date): $1" 1>&2
}
if [[ "$output_vcf" == "" ]] ; then
	write_status "$USAGE_MESSAGE"
	write_status "Output VCF not specified. Use --output to specify output file."
	exit $EX_USAGE
fi
##### --workingdir
write_status "Using working directory \"$workingdir\""
if [[ "$workingdir" == "" ]] ; then
	$working_dir="$(dirname $output_vcf)"
fi
if [[ "$(tr -d ' 	\n' <<< "$workingdir")" != "$workingdir" ]] ; then
		write_status "workingdir cannot contain whitespace"
		exit $EX_USAGE
	fi
workingdir=$(dirname $workingdir/placeholder)
workingdir=$workingdir/$(basename $output_vcf).virusbreakend.working
if [[ ! -d $workingdir ]] ; then
	mkdir -p $workingdir
	if [[ ! -d $workingdir ]] ; then
		write_status "Unable to create $workingdir"
		exit $EX_CANTCREAT
	fi
fi
timestamp=$(date +%Y%m%d_%H%M%S)
# Logging
logfile=$workingdir/full.$timestamp.$HOSTNAME.$$.log
# $1 is message to write
write_status() { # After logging initialised
	echo "$(date): $1" | tee -a $logfile 1>&2
}
write_status "Full log file is: $logfile"
# Timing instrumentation
timinglogfile=$workingdir/timing.$timestamp.$HOSTNAME.$$.log
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
mkdir -p $(dirname $output_vcf)
if [[ ! -d $(dirname $output_vcf) ]] ; then
	write_status "Unable to create directory for $output_vcf for output VCF."
	exit $EX_CANTCREAT
fi
write_status "Using output VCF $output_vcf"
##### --threads
if [[ "$threads" -lt 1 ]] ; then
	write_status "$USAGE_MESSAGE"
	write_status "Illegal thread count: $threads. Specify an integer thread count using the --threads command line argument"
	exit $EX_USAGE
fi
write_status  "Using $threads worker threads."
if [[ "$@" == "" ]] ; then
	write_status  "$USAGE_MESSAGE"
	write_status  "At least one input bam must be specified."
fi
if [[ $force != "true" ]] ; then
	for f in "$@" ; do
		if [[ ! -f $f ]] ; then
			write_status "Input file $f does not exist"
			exit $EX_NOINPUT
		fi
	done
fi
if [[ "$kraken2db" == "" ]] ; then
	echo "$USAGE_MESSAGE"
	write_status "Missing Kraken2 database location. Specify with --kraken2db"
	exit $EX_USAGE
fi
if [[ ! -d "$kraken2db" ]] ; then
	echo "$USAGE_MESSAGE"
	write_status "Unable to find kraken2 database directory '$kraken2db'" 
	exit $EX_NOINPUT
fi
if [[ "$nodesdmp" == "" ]] ; then
	nodesdmp="$kraken2db/taxonomy/nodes.dmp"
fi
if [[ ! -f "$nodesdmp" ]] ; then
	echo "$USAGE_MESSAGE"
	write_status "Unable to find NCBI nodes.dmp file. Specify with --nodesdmp."
	write_status "kraken2-build will include this file in taxonomy/nodes.dmp if --clean was not run."
	write_status "nodes.dmp can be downloaded from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip"
	exit $EX_NOINPUT
fi
taxid_args="-TAXONOMY_IDS null"
if [[ "$hosttaxid" -gt 0 ]] ; then
	if [[ ! -f "$virushostdb" ]] ; then
		write_status "Missing virushostdb file. Specify with --virushostdb"
		write_status "virushostdb can be downloaded from ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv"
	fi
	write_status "Found virushostdb file $virushostdb"
	i=0
	for taxid in $(cut -f 1,8 $virushostdb | grep "	$hosttaxid\$" | cut -f 1) ; do
		i=$(($i + 1))
		taxid_args="$taxid_args -TAXONOMY_IDS $taxid"
	done
	write_status "Found $i viral sequences associated with host NCBI:txid${hosttaxid}"
	if [[ $i -eq 0 ]]; then
		write_status "Terminating early due to 0 sequences of interest for host NCBI:txid${hosttaxid}"
		exit $EX_CONFIG
	fi
else
	taxid_args="$taxid_args -TAXONOMY_IDS 10239" 	# All viruses
fi
if [[ $force != "true" ]] ; then
	for f in "$@" ; do
		if [[ "$(tr -d ' 	\n' <<< "$f")" != "$f" ]] ; then
			write_status "input filenames and paths cannot contain whitespace"
			exit $EX_USAGE
		fi
		write_status "Using input file $f"
	done
fi
# Validate tools exist on path
for tool in kraken2 gridss.sh gridss_annotate_vcf_kraken2.sh samtools java bwa Rscript ; do
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
	write_status "java version: $(java -version 2>&1 | tr '\n' '\t')"
else
	write_status "Unable to run GRIDSS jar. Java 1.8 or later is required."
	write_status "java version: $(java -version  2>&1)"
	exit $EX_CONFIG
fi

# Check kraken2 library files
library_arg=""
for fna in $(find $kraken2db -name library.fna) ; do
	if [[ ! -f $fna.fai ]] ; then
		write_status "Indexing $fna (once-off operation)"
		samtools faidx $fna
	fi
	library_arg="$library_arg -KRAKEN_REFERENCES $fna"
done
if [[ "$library_arg" == "" ]] ; then
	write_status "Unable to find any library.fna files in '$kraken2db'."
	write_status "VIRUSbreakend requires the viral kraken2 reference genomes to be retained."
	write_status "Download using \'kraken2-build --download-library viral --db \"$kraken2db\"'"
	write_status "and do not run kraken2-build --clean as it will remove these files."
	exit $EX_NOINPUT
fi

ulimit -n $(ulimit -Hn) # Reduce likelihood of running out of open file handles 
unset DISPLAY # Prevents errors attempting to connecting to an X server when starting the R plotting device

jvm_args=" \
	-Dpicard.useLegacyParser=false \
	-Dsamjdk.use_async_io_read_samtools=true \
	-Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dsamjdk.buffer_size=4194304 \
	-Dsamjdk.async_io_read_threads=$threads"

file_prefix=$workingdir/$(basename $output_vcf)
file_assembly=$file_prefix.assembly.bam
file_gridss_vcf=$file_prefix.gridss.vcf
file_host_annotated_vcf=$file_prefix.gridss.host_annotated.vcf
file_kraken_annotated_vcf=$file_prefix.gridss.fully_annotated.vcf
file_readname=$file_prefix.readnames.txt
file_report=$file_prefix.kraken2.report.all.txt
file_viral_report=$file_prefix.kraken2.report.viral.txt
file_viral_fa=$file_prefix.viral.fa
exec_concat_fastq=$file_prefix.cat_input_as_fastq.sh

# TODO: Add fastq/fasta version of the pipeline
if [[ ! -f $file_readname ]] ; then
	write_status "Identifying viral sequences"
	rm -f $file_prefix.readnames.txt.tmp
	cat > $exec_concat_fastq << EOF
#!/bin/bash
java -Xmx256m $jvm_args -Dgridss.async.buffersize=2048 -cp $gridss_jar gridss.UnmappedSequencesToFastq \\
EOF
	for f in "$@" ; do
		echo "	-INPUT $f \\" >> $exec_concat_fastq
	done
	cat >> $exec_concat_fastq << EOF
	-OUTPUT /dev/stdout \\
	-INCLUDE_SOFT_CLIPPED_BASES true \\
	-MIN_SEQUENCE_LENGTH 20 \\
	-UNIQUE_NAME false
EOF
	chmod +x $exec_concat_fastq
	{ $timecmd $exec_concat_fastq \
	| kraken2 \
		--threads $threads \
		--db $kraken2db \
		--report $file_report \
		$kraken2args \
		/dev/stdin \
	| java -Xmx512m $jvm_args -cp $gridss_jar gridss.kraken.SubsetToTaxonomy \
		-INPUT /dev/stdin \
		-OUTPUT $file_readname.tmp \
		-FORMAT READ_NAME \
		-NCBI_NODES_DMP $nodesdmp \
	&& mv $file_readname.tmp $file_readname \
	; } 1>&2 2>> $logfile
fi
if [[ ! -f $file_viral_fa ]] ; then
	write_status "Identifying viruses in sample"
	{ $timecmd java -Xmx4g $jvm_args -cp $gridss_jar gridss.kraken.ExtractBestSequencesBasedOnReport \
		-INPUT $file_report \
		-OUTPUT $file_viral_fa \
		-REPORT_OUTPUT $file_viral_report \
		-NCBI_NODES_DMP $nodesdmp \
		-KRAKEN_REFERENCES $kraken2db/library/viral/library.fna \
		-MIN_SUPPORTING_READS $minreads \
		-SEQUENCES_TO_RETURN $viralgenomes \
		$taxid_args \
	; } 1>&2 2>> $logfile
fi
if [[ ! -s $file_viral_fa ]] ; then
	write_status "No viral sequences supported by at least $minreads reads."
	trap - EXIT
	exit 0 # success!
fi
if [[ ! -f $file_viral_fa.bwt ]] ; then
	write_status "Creating index of viral sequences"
	{ $timecmd samtools faidx $file_viral_fa && bwa index $file_viral_fa ; } 1>&2 2>> $logfile
fi

gridss_input_args=""
for f in "$@" ; do
	basef=$(basename "$f" | tr -d ' 	<(:/|$@&%^\\)>')
	infile_prefix=$workingdir/$basef
	infile_fq=$infile_prefix.viral.unpaired.fq
	infile_fq1=$infile_prefix.viral.R1.fq
	infile_fq2=$infile_prefix.viral.R2.fq
	infile_unsorted_sam=$infile_prefix.viral.sam
	infile_bam=$infile_prefix.viral.bam
	
	if [[ ! -f $infile_fq ]] ; then
		exec_concat_fastq=$infile_prefix.extract_reads.sh
		write_status "Extracting viral reads	$f"
		cat > $exec_concat_fastq << EOF
java -Xmx4g $jvm_args -cp $gridss_jar gridss.ExtractFragmentsToFastq \\
	-INPUT $f \\
	-READ_NAMES $file_readname \\
	-OUTPUT_FQ $infile_fq \\
	-OUTPUT_FQ1 $infile_fq1 \\
	-OUTPUT_FQ2 $infile_fq2
EOF
		chmod +x $exec_concat_fastq
		{ $timecmd $exec_concat_fastq ; } 1>&2 2>> $logfile
	fi
	if [[ ! -f $infile_bam ]] ; then
		write_status "Aligning viral reads	$f"
		{ $timecmd cat \
			<(bwa mem -Y -t $threads $file_viral_fa $infile_fq1 $infile_fq1) \
			<(bwa mem -Y -t $threads $file_viral_fa $infile_fq | grep -v "^@") \
		| samtools fixmate -m -O BAM - - \
		| samtools sort -l 0 -@ $threads -T $infile_unsorted_sam.tmp.sorted - \
		| samtools markdup -O BAM -@ $threads - $infile_bam.tmp.bam \
		&& mv $infile_bam.tmp.bam $infile_bam \
		&& samtools index $infile_bam \
		; } 1>&2 2>> $logfile
	fi
	gridss_dir=$workingdir/$(basename $infile_bam).gridss.working
	gridss_prefix=gridss_dir/$(basename $infile_bam)
	if [[ ! -f $gridss_prefix.insert_size_metrics ]] ; then
		write_status "Gathering metrics from host alignment	$f"
		# Ideally the metrics on the viral sequence would match the metrics from the host.
		# Unfortunately, this is generally not the case.
		# To ensure we assemble fragments correctly, we need to grab the host alignment metrics.
		# If GRIDSS has been run, we could use that but we don't want GRIDSS to be an explicit
		# requirement of this pipeline
		# This approach doesn't work for fastq input files.
		mkdir -p $gridss_dir
		{ $timecmd java -Xmx4g $jvm_args \
			-cp $gridss_jar gridss.analysis.CollectGridssMetrics \
			--INPUT $infile_bam \
			--OUTPUT $gridss_prefix \
			--REFERENCE_SEQUENCE $file_viral_fa \
			--THRESHOLD_COVERAGE $metricsmaxcoverage \
			--TMP_DIR $workingdir \
			--FILE_EXTENSION null \
			--STOP_AFTER $metricsrecords \
		; } 1>&2 2>> $logfile
	fi
	gridss_input_args="$gridss_input_args $infile_bam"
done
if [[ ! -f $file_gridss_vcf ]] ; then
	write_status "Calling structural variants"
	{ $timecmd gridss.sh \
		-w $workingdir \
		-t $threads \
		-r $file_viral_fa \
		-j $gridss_jar \
		-o $file_gridss_vcf \
		-a $file_assembly \
		--maxcoverage $maxcoverage \
		$gridssargs \
		$gridss_input_args \
	; } 1>&2 2>> $logfile
fi
if [[ ! -f $file_host_annotated_vcf ]] ; then
	write_status "Annotating host genome integrations"
	{ $timecmd java -Xmx4g $jvm_args \
			-Dgridss.output_to_temp_file=true \
			-cp $gridss_jar gridss.AnnotateInsertedSequence \
			--TMP_DIR $workingdir \
			--WORKING_DIR $workingdir \
			--REFERENCE_SEQUENCE $reference \
			--WORKER_THREADS $threads \
			--INPUT $file_gridss_vcf \
			--OUTPUT $file_host_annotated_vcf \
	; } 1>&2 2>> $logfile
fi
if [[ ! -f $file_kraken_annotated_vcf ]] ; then
	write_status "Annotating kraken2 host"
	{ $timecmd gridss_annotate_vcf_kraken2.sh \
		-o $file_kraken_annotated_vcf \
		-j $gridss_jar \
		--kraken2db $kraken2db \
		--threads $threads \
		$kraken2args \
		$file_host_annotated_vcf  \
	; } 1>&2 2>> $logfile
fi

cp $file_kraken_annotated_vcf $output_vcf

write_status "Generated $output_vcf"
write_status "Done"

trap - EXIT
exit 0 # success!
