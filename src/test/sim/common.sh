#!/bin/bash
#
# Common variables and utility functions
#
CONFIG=fastcompare
#CONFIG=fullmatrix

if [[ -d data.$1 ]] ; then
	echo "using data.$1"
	CONFIG=$1
fi

EXECUTION_CONTEXT=torque #can be one of {local, torque, slurm}
if [[ -f "/usr/local/bioinf/newpaths" ]] ; then
	# use Java 1.8
	source /usr/local/bioinf/newpaths
fi

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#HARNESS_CLASSPATH="$CLASSPATH:$BASE_DIR/tools/harness:$(find $BASE_DIR/tools/picard $BASE_DIR/tools/harness -name '*.jar' | tr '\n' ':')"
HARNESS_CLASSPATH=~/dev/ws/ViralInsertionDetection/bin:~/dev/ws/ViralInsertionDetection/dist/*:~/dev/ws/ViralInsertionDetection/lib/*:~/dev/ws/ViralInsertionDetection/lib/biojava/*:~/dev/ws/ViralInsertionDetection/lib/picard/*:$CLASSPATH
DATA_DIR=$BASE_DIR/data.$CONFIG
mkdir -p $DATA_DIR

STARTING_DEPTH=100
READ_DEPTHS="$STARTING_DEPTH 60 30 15 8 4"
READ_LENGTHS="36 50 75 100 150 250"
FRAGMENT_SIZE="500 400 300 250 200 150"
FULL_MATRIX=0

MIN_EXCLUDED_CALL_DEPTH=130 #temp hack so we're not waiting 3 weeks for full 100x runs of everything
# override default settings
source settings.$CONFIG

trap "exit 1" TERM
export TOP_PID=$$
function ls_reference_vcf {
	for FILE in $DATA_DIR/*.reference.vcf ; do
		echo $FILE
	done
}
# $1 sort order: sq for queryname sorted, sc for coordinate sorted
# $2 aligner: 
function ls_aligned_bam {
	for FILE in $DATA_DIR/*.${1:-sc}.bam; do
		if [ "$2" != "" ] ; then
			if grep -q -E "^CX_ALIGNER=$2\$" $(get_cx_prefix $FILE).metadata; then
				true # we're good
			else
				continue
			fi
		fi
		echo $FILE
	done
}
# includes all context variables CX_*
function cx_save {
	local METADATA_TMP_FILE=$(mktemp)
	set | grep -E '^CX_' > $METADATA_TMP_FILE
	CX=$DATA_DIR/$(md5sum $METADATA_TMP_FILE | cut -b -32)
	if [ ! -s $CX.metadata ] ; then
		cp $METADATA_TMP_FILE $CX.metadata
	fi
}
function try_get_cx_prefix {
	local NAME=$1
	if [[ $NAME =~ .*/[0-9a-f]{32} ]] ; then
		NAME=${BASH_REMATCH[0]}
		while [[ ${NAME%/*} =~ .*/[0-9a-f]{32} ]] ; do
			NAME=${NAME%/*}
		done
		echo $NAME
		return
	fi
	echo ""
}
function get_cx_prefix {
	local PREFIX=$(try_get_cx_prefix $1)
	if [ ! -f "$PREFIX.metadata" ] ; then
		echo "Error: missing metadata file for $1" 1>&2
		kill -s TERM $TOP_PID; exit 1
	fi
	if [ "$PREFIX" == "" ] ; then
		echo "Prefix not found for $1" 1>&2
		kill -s TERM $TOP_PID; exit 1
	fi
	echo $PREFIX
}

# changes the metadata context to the given context
function cx_load {
	PREFIX=$(get_cx_prefix $1)
	cx_flush
	if [ "$PREFIX" == "" ] ; then
		echo "Prefix not found for $1" 1>&2
		kill -s TERM $TOP_PID; exit 1
	fi
	# load context from the metadata file
	. $PREFIX.metadata
	CX="$PREFIX"
	if [[ "$CX_READ_FRAGMENT_LENGTH" != "" ]] ; then
		if [[ "$CX_READ_FRAGMENT_STDDEV" == "" ]] ; then
			CX_READ_FRAGMENT_STDDEV=$(( CX_READ_FRAGMENT_LENGTH / 10))
		fi
	fi
}
function cx_flush {
	# removes all variables starting with CX_
	unset -v $(set | grep -E '^CX_' | cut -s -d '=' -f 1)
}
function xc_flush {
	# removes all variables starting with CX_
	unset -v $(set | grep -E '^XC_' | cut -s -d '=' -f 1)
	# requires cluster head to match CPU specs of cluster nodes
	XC_CORES=$MAX_CORES
	XC_MEMORY=$MAX_MEMORY
}
# XC - execution context
# when running any given script, instead of running it directly
# from the shell script, XC_SCRIPT should be populated with
# the commands to execute.
# XC Variants
# XC_SCRIPT commands to execute
# XC_TRAP clean-up commands to execute when exiting
# XC_CORES number of core available
# XC_MULTICORE set iff process is multithreaded
# XC_NOCLUSTER set iff process is unable to run on the cluster
# XC_MEMORY expected maximum memory in MB
# XC_OUTPUT output file. If this file exists, the script will not be rerun
MAX_CORES=$(grep processor < /proc/cpuinfo | wc -l)
MAX_MEMORY=$(free -m | gawk '{ if (NR == 2) { print $2 } } ')

# Executes the current execution script
function xc_exec {
	if [ "$CX" == "" ] ; then
		echo "Error: CX not set when calling xc_exec" 1>&2
		kill -s TERM $TOP_PID; exit 1
	fi
	if [ ! -s "$XC_OUTPUT" ] ; then
		if [ -f $CX.lock ]  ; then
			echo "$(basename $CX): Not scheduling as lock file exists ($(cat $CX.lock))"
		#TODO: check for any upstream locks based on prefix match to CX_*
		elif [[ "$CX_CALLER" != "" && "$CX_READ_DEPTH" -ge $MIN_EXCLUDED_CALL_DEPTH ]] ; then
			echo "$(basename $CX): Not scheduling ${CX_READ_DEPTH}x coverage as it takes too long"
		else
			echo "$(tput setaf 2)$(basename $CX): Scheduling task$(tput sgr0)"
			xc_exec_$EXECUTION_CONTEXT
		fi
	else
		echo "$(basename $CX): found $XC_OUTPUT. Task does not need to be rerun"
	fi
	xc_flush
}
# $1 expected lock contents
function write_xc_exec_scripts {
	cat > $CX$XC_SUFFIX.script.sh << EOF
#!/bin/bash
$XC_SCRIPT
EOF
	local RES="mem=${XC_MEMORY}mb"
	if [[ "$XC_MULTICORE" != "" || $((XC_MEMORY * 2)) -gt $MAX_MEMORY ]] ; then
		# reserve the entire node for multicore or high-memory jobs
		RES="nodes=1:ppn=${MAX_CORES}"
	fi
	# TODO: XC_CORES and XC_MEMORY
	cat > $CX$XC_SUFFIX.sh << EOF
#!/bin/bash
#PBS -S /bin/bash
#PBS -e $CX$XC_SUFFIX.stderr
#PBS -o $CX$XC_SUFFIX.stdout
#PBS -N indel$(basename $CX)
#PBS -V
#PBS -l $RES
set -o pipefail
#renice 20 -p \$\$

if [[ "\$(cat $CX.lock)" != "$LOCKKEY" ]] ; then
		echo "Found unexpected lock key \"\$(cat $CX.lock)\" for $(basename $CX): not processing"
		exit 1
fi
if pidof $CX$XC_SUFFIX.script.sh ; then
	echo "Found running process $CX$XC_SUFFIX.script.sh: not processing"
	exit 1
fi
killall -p $CX$XC_SUFFIX.script.sh 2>/dev/null
trap "{ rm -f $CX.lock ; ${XC_TRAP:-true}; exit 255; }" EXIT
cd "$BASE_DIR"
/usr/bin/time -f "user	%U	elapsed	%e	sys	%S	maxResMemKb	%M	avgTotMemKb	%K I	%I	O	%O	exit	%x" -o $CX$XC_SUFFIX.time $CX$XC_SUFFIX.script.sh
EOF
	chmod +x $CX$XC_SUFFIX.sh $CX$XC_SUFFIX.script.sh
}
# Executes the current execution script locally
function xc_exec_local {
	if [[ -f $CX.lock ]] ; then
		echo "Found unexpected lock for $(basename $CX): not processing"
		return
	fi
	LOCKKEY="$(hostname):$TOP_PID"
	echo $LOCKKEY > $CX.lock
	write_xc_exec_scripts
	chmod a+x $CX$XC_SUFFIX.sh
	( $CX$XC_SUFFIX.sh >(tee $CX$XC_SUFFIX.stdout) 2> >(tee $CX$XC_SUFFIX.stderr >&2) )
	# run the clean up code since trap is ineffective when we source the script
	rm -f $CX.lock; ${XC_TRAP:-true}
}
# Executes the current execution script on a cluster
function xc_exec_torque {
	if $(which qsub > /dev/null 2>&1) ; then
		if [ "$XC_NOCLUSTER" != "" ] ; then
			echo "$(tput setaf 3)$(basename $CX): Not scheduling as it has been flagged NOCLUSTER$(tput sgr0)"
			xc_flush
			return
		fi
		# use exclusive nodes for now as we haven't profiled memory usage
		#if [ "$XC_MULTICORE" != "" ] ; then
		write_xc_exec_scripts
		rm -f $CX$XC_SUFFIX.stderr $CX$XC_SUFFIX.stdout
		qsub $CX$XC_SUFFIX.sh | tee $CX.lock
	else
		# fall back to local execution if qsub not available
		echo "qsub not found. Executing locally."
		xc_exec_local
	fi
}

cx_flush
xc_flush

