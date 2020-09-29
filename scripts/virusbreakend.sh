#!/bin/bash
# ../scripts/virusbreakend.sh -j ../target/gridss-2.10.0-gridss-jar-with-dependencies.jar -o vbe_out.vcf -r ../../ref/hg19.fa --db ../../virusbreakend/virusbreakenddb ERR093636_virusbreakend_minimal_example.bam
# ../scripts/virusbreakend.sh -j ../target/gridss-2.10.0-gridss-jar-with-dependencies.jar -o vbe_out.vcf -r ../../ref/hg19.fa --db ../../virusbreakend/virusbreakenddb ERR093636_virusbreakend_minimal_example_slower_fastq_input_R1.fq ERR093636_virusbreakend_minimal_example_slower_fastq_input_R2.fq
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
metricsrecords=1000000
metricsmaxcoverage=100000
maxcoverage=1000000
hosttaxid=9606
force="false"
forceunpairedfastq="false"
USAGE_MESSAGE="
VIRUSBreakend: Viral Integration Recognition Using Single Breakends

Usage: virusbreakend.sh [options] input.bam

	-r/--reference: reference genome of host species.
	-o/--output: output VCF.
	-j/--jar: location of GRIDSS jar
	-t/--threads: number of threads to use. (Default: $threads).
	-w/--workingdir: directory to place intermediate and temporary files. (Default: $workingdir).
	--db: path to virusbreakenddb database directory. Use the supplied virusbreakend-build.sh to build.
	--hosttaxid: NCBI taxonomy id of host. Used to filter viral sequences of interest to those infecting this host. Default: $hosttaxid)
	--kraken2args: additional kraken2 arguments
	--gridssargs: additional GRIDSS arguments
	--minreads: minimum number of viral reads perform integration detection (Default: $minreads)
	--viralgenomes: number of viral genomes to consider. Multiple closely related genomes will result in a high false negative rate due to multi-mapping reads. (Default: $viralgenomes)
	"
# handled by virusbreakend-build.sh
#--kraken2db: kraken2 database
#--virushostdb: location of virushostdb.tsv. Available from ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv (Default: {kraken2db}/virushostdb.tsv)
#--nodesdmp: location of NCBI nodes.dmp. Can be downloaded from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip. (Default: {kraken2db}/taxonomy/nodes.dmp)
OPTIONS=ho:t:j:w:r:f
LONGOPTS=help,output:,jar:,threads:,reference:,workingdir:,db:,kraken2db:,kraken2args:,gridssargs:,nodesdmp:,minreads:,hosttaxid:,virushostdb:,force,forceunpairedfastq
! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
    #  then getopt has complained about wrong arguments to stdout
	echo "$USAGE_MESSAGE" 1>&2
    exit $EX_USAGE
fi
eval set -- "$PARSED"
while true; do
	case "$1" in
		-h|--help)
			echo "$USAGE_MESSAGE" 1>&2
			exit 0
			;;
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
		--db|--kraken2db)
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
		--forceunpairedfastq)
			forceunpairedfastq="true"
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
logfile=$workingdir/virusbreakend.$timestamp.$HOSTNAME.$$.log
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
	exit $EX_USAGE
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
if [[ "$virushostdb" == "" ]] ; then
	virushostdb="$kraken2db/virushostdb.tsv"
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
for tool in kraken2 gridss.sh gridss.sh gridss_annotate_vcf_kraken2.sh gridss_annotate_vcf_repeatmasker.sh samtools java bwa Rscript ; do
	if ! which $tool >/dev/null; then
		write_status "Error: unable to find $tool on \$PATH"
		exit $EX_CONFIG
	fi
	write_status "Found $(which $tool)"
done
if which gridsstools > /dev/null ; then
	write_status "Found $(which gridsstools)"
	write_status "gridsstools version: $(gridsstools --version)"
else 
	write_status "MISSING gridsstools. Execution will take 2-3x time longer than when using gridsstools."
fi
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


# Hack to support streaming using bash redirection
# "<(cat file.fastq)" notation for fastq
# $1: filename or bash redirection string
function clean_filename {
	# handle file descriptor redirection
	basename "$1" | tr -d ' 	<(:/|$@&%^\\)>'
}

function echo_cat_uncompressed {
	cleanf=$(clean_filename "$1")
	if [[ "${cleanf%.gz}" != "$cleanf" ]] ; then
		if which pigz > /dev/null ; then
			echo -n "pigz -c -d -p $threads $1"
		else
			write_status "WARNING: pigz not found on PATH. Using single-threaded gunzip for $1"
			echo -n "gunzip -c $1"
		fi
	else
		echo -n "cat $1"
	fi
}

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
file_host_annotated_vcf=$file_prefix.gridss.host.vcf
file_kraken_annotated_vcf=$file_prefix.gridss.host.k2.vcf
file_rm_annotated_vcf=$file_prefix.gridss.host.k2.rm.vcf
file_filtered_vcf=$file_prefix.gridss.host.k2.rm.filtered.vcf
file_readname=$file_prefix.readnames.txt
file_report=$file_prefix.kraken2.report.all.txt
file_viral_report=$file_prefix.kraken2.report.viral.txt
file_extracted_report=$file_prefix.kraken2.report.viral.extracted.txt
file_viral_fa=$file_prefix.viral.fa
file_merged_bam=$file_prefix.merged.bam
file_wgs_metrics=$file_prefix.wgs_metrics.txt
exec_concat_fastq=$file_prefix.cat_input_as_fastq.sh

if [[ ! -f $file_readname ]] ; then
	write_status "Identifying viral sequences"
	rm -f $exec_concat_fastq $file_prefix.readnames.txt.tmp
	echo "#!/bin/bash" > $exec_concat_fastq
	for f in "$@" ; do
		cleanf=$(clean_filename "$f")
		if [[ "${cleanf%.fastq.gz}" != "$cleanf" ]] ||
				[[ "${cleanf%.fq.gz}" != "$cleanf" ]] ||
				[[ "${cleanf%.fq}" != "$cleanf" ]] || 
				[[ "${cleanf%.fastq}" != "$cleanf" ]] ; then
			# fastq input
			echo_cat_uncompressed "$f" >> $exec_concat_fastq
			echo >> $exec_concat_fastq
		else
			# BAM input
			if which gridsstools > /dev/null ; then
				echo "gridsstools unmappedSequencesToFastq -@ $threads $f" >> $exec_concat_fastq
			else
				cat >> $exec_concat_fastq << EOF
java -Xmx256m $jvm_args -Dgridss.async.buffersize=2048 -cp $gridss_jar gridss.UnmappedSequencesToFastq \\
	-INPUT $f \\
	-OUTPUT /dev/stdout \\
	-INCLUDE_SOFT_CLIPPED_BASES true \\
	-MIN_SEQUENCE_LENGTH 20 \\
	-UNIQUE_NAME false
EOF
			fi
		fi
	done
	chmod +x $exec_concat_fastq
	{ $timecmd $exec_concat_fastq \
	| kraken2 \
		--threads $threads \
		--db $kraken2db \
		--report $file_report \
		$kraken2args \
		/dev/stdin \
	| java -Xmx512m $jvm_args -cp $gridss_jar gridss.kraken.SubsetToTaxonomy \
		--INPUT /dev/stdin \
		--OUTPUT $file_readname.tmp \
		--FORMAT READ_NAME \
		--NCBI_NODES_DMP $nodesdmp \
	&& mv $file_readname.tmp $file_readname \
	; } 1>&2 2>> $logfile
else
	write_status "Identifying viral sequences	Skipped: found	$file_readname"
fi
if [[ ! -f $file_viral_fa ]] ; then
	write_status "Identifying viruses in sample based on kraken2 summary report"
	# The sort is so we will include the virushostdb sequences in
	# library/added before the kraken sequences in library/viral
	kraken_references_arg=$(for fa in $(find $kraken2db -path '**/library/**/*.fna' | sort) ; do echo -n "--KRAKEN_REFERENCES $fa "; done)
	{ $timecmd java -Xmx4g $jvm_args -cp $gridss_jar gridss.kraken.ExtractBestSequencesBasedOnReport \
		--INPUT $file_report \
		--OUTPUT $file_viral_fa \
		--REPORT_OUTPUT $file_viral_report \
		--SUMMARY_REPORT_OUTPUT $file_extracted_report \
		--NCBI_NODES_DMP $nodesdmp \
		$kraken_references_arg \
		--KRAKEN_REFERENCES $kraken2db/library/viral/library.fna \
		--MIN_SUPPORTING_READS $minreads \
		--TAXID_TO_RETURN $viralgenomes \
		$taxid_args \
	; } 1>&2 2>> $logfile
else
	write_status "Identifying viruses	Skipped: found	$file_viral_fa"
fi
if [[ ! -s $file_viral_fa ]] ; then
	write_status "No viral sequences supported by at least $minreads reads."
	trap - EXIT
	exit 0 # success!
fi
if [[ ! -f $file_viral_fa.bwt ]] ; then
	write_status "Creating index of viral sequences"
	{ $timecmd samtools faidx $file_viral_fa && bwa index $file_viral_fa ; } 1>&2 2>> $logfile
else
	write_status "Creating index of viral sequences	Skipped: found	$file_viral_fa.bwt"
fi

bam_list_args=""
for f in "$@" ; do
	cleanf=$(clean_filename "$f")
	infile_prefix=$workingdir/$cleanf
	infile_fq=$infile_prefix.viral.unpaired.fq
	infile_fq1=$infile_prefix.viral.R1.fq
	infile_fq2=$infile_prefix.viral.R2.fq
	infile_unsorted_sam=$infile_prefix.viral.sam
	infile_bam=$infile_prefix.viral.bam
	fastq_extension=""
	for possible_extension in .fastq.gz .fq.gz .fq .fastq ; do
		if [[ $(basename $cleanf $possible_extension) != "$cleanf" ]] ; then
			fastq_extension=$possible_extension
		fi
	done
	fq2=""
	fq1=""
	if [[ "$fastq_extension" != "" ]] ; then
		# fastq handling
		for f2 in "$@" ; do
			cleanf2=$(clean_filename "$f2")
			if [[ $(basename "$cleanf2" 2$fastq_extension)1$fastq_extension == $cleanf ]] ; then
				fq1=$f
				fq2=$f2
				write_status "Treating as fastq pair	$f	$f2"
			fi
			if [[ $(basename "$cleanf2" 1$fastq_extension)2$fastq_extension == $cleanf ]] ; then
				fq1=$f2
				fq2=$f
				write_status "Treating as second of fastq pair	$f"
			fi
		done
		if [[ "$f" == "$fq2" ]] ; then
			# this is f2 in a paired fastq - skip processing since we do it with the first in pair
			continue
		fi
		if [[ "$fq1" == "" ]] ; then
			write_status "Could not find a pairing for $1 based on filename suffix."
			if [[ "$forceunpairedfastq" != true ]] ; then
				write_status "VIRUSbreakend expects paired fastq inputs to have \"1\" and \"2\" immediately before the fastq suffix."
				write_status "For example, input_1.fastq.gz, input_2.fastq.gz or input_R1.fq, input_R2.fq"
				write_status "If this file is really single-end short read sequencing data, use --forceunpairedfastq to process as unpaired"
				exit EX_CONFIG
			fi
		fi
	fi
	if [[ ! -f $infile_fq ]] ; then
		exec_extract_reads=$infile_prefix.extract_reads.sh
		write_status "Extracting viral reads	$f	$fq2"
		rm -f $exec_extract_reads
		echo "#!/bin/bash" > $exec_extract_reads
		if [[ "$fastq_extension" == "" ]] ; then
			if which gridsstools > /dev/null ; then
				cat >> $exec_extract_reads << EOF
gridsstools extractFragmentsToFastq \
	-@ $threads \
	-r $file_readname \
	-o $infile_fq \
	-1 $infile_fq1 \
	-2 $infile_fq2 \
	$f
EOF
			else
				cat >> $exec_extract_reads << EOF
java -Xmx4g $jvm_args -cp $gridss_jar gridss.ExtractFragmentsToFastq \\
	-INPUT $f \\
	-READ_NAMES $file_readname \\
	-OUTPUT_FQ $infile_fq \\
	-OUTPUT_FQ1 $infile_fq1 \\
	-OUTPUT_FQ2 $infile_fq2
EOF
			fi
		else
			# TODO: replace awk with C code: https://www.biostars.org/p/10353/
			awk_script='NR==FNR { lookup["@" $1]=1; next } { if (lookup[$1] == 1) { print; getline; print; getline; print; getline; print } else { getline ; getline; getline } }'
			if [[ "$fq1" != "" ]] ; then
				# AWK script to pull out
				echo "$(echo_cat_uncompressed $fq1) | awk '$awk_script' $file_readname /dev/stdin > $infile_fq1 &" >> $exec_extract_reads
				echo "$(echo_cat_uncompressed $fq2) | awk '$awk_script' $file_readname /dev/stdin > $infile_fq2 &" >> $exec_extract_reads
				echo "echo -n > $infile_fq &" >> $exec_extract_reads
				echo "wait" >> $exec_extract_reads
			else
				echo "$(echo_cat_uncompressed $f) | awk '$awk_script' $file_readname /dev/stdin > $infile_fq &" >> $exec_extract_reads
				echo "echo -n > $infile_fq1 &" >> $exec_extract_reads
				echo "echo -n > $infile_fq2 &" >> $exec_extract_reads
				echo "wait &" >> $exec_extract_reads
			fi
		fi
		chmod +x $exec_extract_reads
		{ $timecmd $exec_extract_reads ; } 1>&2 2>> $logfile
	else
		write_status "Extracting viral reads	Skipped: found	$infile_fq"
	fi
	if [[ ! -f $infile_bam ]] ; then
		write_status "Aligning viral reads	$f"
		{ $timecmd cat \
			<(bwa mem -Y -t $threads $file_viral_fa $infile_fq1 $infile_fq2) \
			<(bwa mem -Y -t $threads $file_viral_fa $infile_fq | grep -v "^@") \
		| samtools fixmate -m -O BAM - - \
		| samtools sort -l 0 -@ $threads -T $infile_unsorted_sam.tmp.sorted - \
		| samtools markdup -O BAM -@ $threads - $infile_bam.tmp.bam \
		&& mv $infile_bam.tmp.bam $infile_bam \
		&& samtools index $infile_bam \
		; } 1>&2 2>> $logfile
	else
		write_status "Aligning viral reads	Skipped: found	$infile_bam"
	fi
	gridss_dir=$workingdir/$(basename $infile_bam).gridss.working
	gridss_prefix=$gridss_dir/$(basename $infile_bam)
	if [[ $fastq_extension == "" ]] ; then
		if [[ ! -f $gridss_prefix.insert_size_metrics ]] ; then
			write_status "Gathering metrics from host alignment	$f"
			# Ideally the metrics on the viral sequence would match the metrics from the host.
			# Unfortunately, this is generally not the case.
			# To ensure we assemble fragments correctly, we need to grab the host alignment metrics.
			# If GRIDSS has been run, we could use that but we don't want GRIDSS to be an explicit
			# requirement of this pipeline
			# This approach doesn't work for fastq input files.
			exec_extract_host_metrics=$infile_prefix.extract_host_metrics.sh
			rm -f $exec_extract_host_metrics
			cat > $exec_extract_host_metrics << EOF
java -Xmx4g $jvm_args \
	-cp $gridss_jar gridss.analysis.CollectGridssMetrics \
	--INPUT $f \
	--OUTPUT $gridss_prefix \
	--REFERENCE_SEQUENCE $reference \
	--THRESHOLD_COVERAGE $metricsmaxcoverage \
	--TMP_DIR $workingdir \
	--FILE_EXTENSION null \
	--STOP_AFTER $metricsrecords
EOF
			chmod +x $exec_extract_host_metrics
			mkdir -p $gridss_dir
			{ $timecmd $exec_extract_host_metrics; } 1>&2 2>> $logfile
		else
			write_status "Gathering metrics from host alignment	Skipped: found	$gridss_prefix.insert_size_metrics"
		fi
	else
		write_status "Unable to use host metrics for fastq $f	Skipping metrics pre-calculation"
	fi
	bam_list_args="$bam_list_args $infile_bam"
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
		$bam_list_args \
	; } 1>&2 2>> $logfile
else
	write_status "Calling structural variants	Skipped: found	$file_gridss_vcf"
fi
if [[ ! -f $file_host_annotated_vcf ]] ; then
	# Make sure we have the appropriate indexes for the host reference genome
	{ $timecmd gridss.sh \
		-r $reference \
		-t $threads \
		-j $gridss_jar \
		-s setupreference \
		-a $file_assembly \
		-o placeholder.vcf \
		$bam_list_args \
	; } 1>&2 2>> $logfile
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
	# external bwa process
	# --ALIGNER_COMMAND_LINE null \
	# --ALIGNER_COMMAND_LINE bwa \
	# --ALIGNER_COMMAND_LINE mem \
	# --ALIGNER_COMMAND_LINE -L 0,0 \
	# --ALIGNER_COMMAND_LINE -t \
	# --ALIGNER_COMMAND_LINE '%3$d' \
	# --ALIGNER_COMMAND_LINE '%2$s' \
	# --ALIGNER_COMMAND_LINE '%1$s' \
else
	write_status "Annotating host genome integrations	Skipped: found	$file_host_annotated_vcf"
fi
if [[ ! -f $file_kraken_annotated_vcf ]] ; then
	write_status "Annotating kraken2"
	{ $timecmd gridss_annotate_vcf_kraken2.sh \
		-o $file_kraken_annotated_vcf \
		-j $gridss_jar \
		--kraken2db $kraken2db \
		--threads $threads \
		$kraken2args \
		$file_host_annotated_vcf  \
	; } 1>&2 2>> $logfile
else
	write_status "Annotating kraken2	Skipped: found	$file_kraken_annotated_vcf"
fi
if [[ ! -f $file_rm_annotated_vcf ]] ; then
	write_status "Annotating RepeatMasker"
	{ $timecmd gridss_annotate_vcf_repeatmasker.sh \
		-w $workingdir \
		-o $file_rm_annotated_vcf \
		-j $gridss_jar \
		--threads $threads \
		$kraken2args \
		$file_kraken_annotated_vcf \
	; } 1>&2 2>> $logfile
else
	write_status "Annotating RepeatMasker	Skipped: found	$file_rm_annotated_vcf"
fi
if [[ ! -f $file_filtered_vcf ]] ; then
	write_status "Filtering to host integrations"
	{ $timecmd java -Xmx64m $jvm_args -cp $gridss_jar gridss.VirusBreakendFilter \
		--INPUT $file_rm_annotated_vcf \
		--OUTPUT $file_filtered_vcf \
		--REFERENCE_SEQUENCE $reference \
	; } 1>&2 2>> $logfile
fi
if [[ ! -f $file_wgs_metrics ]] ; then
	write_status "Calculating virus WGS metrics"
	{ $timecmd samtools merge -@ $threads $file_merged_bam $bam_list_args && \
		java -Xmx1g $jvm_args \
			-cp $gridss_jar picard.cmdline.PicardCommandLine CollectWgsMetrics \
			--INPUT $file_merged_bam \
			--OUTPUT $file_wgs_metrics \
			--REFERENCE_SEQUENCE $file_viral_fa \
			--COVERAGE_CAP 10000 \
			--COUNT_UNPAIRED true \
	; } 1>&2 2>> $logfile
fi
cp $file_filtered_vcf $output_vcf
cp $file_extracted_report $output_vcf.kraken2.summary.csv
cp $file_wgs_metrics $output_vcf.wgs_metrics.txt

write_status "Generated $output_vcf"
write_status "Done"

trap - EXIT
exit 0 # success!
