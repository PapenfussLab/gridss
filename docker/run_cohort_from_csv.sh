#!/bin/bash
#
# GRIDSS docker pipeline that runs a cohort based on the files listed in a CSV
#
# Lines starting with # are ignored
#
# The columns in the CSV are as follows:
# 1) patient id
# 2) GRIDSS output SV vcf file. A default filename of PATIENT_ID.gridss.vcf will be used if this is blank
# 3) GRIDSS output assembly bam. A default filename of PATIENT_ID.assembly.bam will be used if this is blank
# 4) first sample BAM from patient (for somatic calling this is the matched normal bam)
# 5+) additional sample BAMs (for somatic calling the primary is typically the first somatic sample)
#
##### CHANGE PATH TO YOUR REFERENCE GENOME
REFERENCE=../../ref/hg19.fa
#####
BLACKLIST=ENCFF001TDO.bed
CONTAINER=gridss/gridss

if [[ ! -f "$1" ]] ; then
	echo "Usage run_cohort_from_csv.sh cohort.csv output_directory"  1>&2
	echo "Missing cohort CSV file" 1>&2
	exit 1
fi

if [[ "$2" == "" ]] ; then
	echo "Usage run_cohort_from_csv.sh cohort.csv output_directory"  1>&2
	echo "Missing output directory" 1>&2
	exit 1
fi
OUTDIR="$2"
mkdir -p "$OUTDIR"
if [[ ! -d "$OUTDIR" ]] ; then
	echo "$OUTDIR is not a directory" 1>&2
	exit 1
fi
if [[ "$BLACKLIST" == "ENCFF001TDO.bed" ]] ; then
	BLACKLIST=$OUTDIR/ENCFF001TDO.bed
	if [[ ! -f $BLACKLIST ]] ; then
		cd $OUTDIR
		wget https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz
		gunzip ENCFF001TDO.bed.gz
		cd -
	fi
fi
# Sanity check reference exists
if [[ ! -f "$REFERENCE" ]] ; then
	echo "Missing reference genome $REFERENCE. Update the REFERENCE variable in the shell script to your hg19 location" 1>&2
	echo "For the example file chr12.1527326.DEL1024.bam, ReorderSam can be used to match to your version of hg19. In the case of this example, only \"chr12\" is required to exist, and difference in alternate contigs can be ignored (using ALLOW_INCOMPLETE_DICT_CONCORDANCE=true)." 1>&2
	echo "For real data, please ensure that all BAM files are aligned to the same reference, and the reference supplied to GRIDSS matches that used for alignment." 1>&2
	exit 1
fi
if [[ ! -f $BLACKLIST ]] ; then
	echo "Missing blacklist file $BLACKLIST." 1>&2
	echo "If you do not wish to use a blacklist, either update this script or supply an empty BED file as the blacklist." 1>&2
	exit 1
fi

# Check that we have a bwa index.
# Technically this check could be done inside a script in the docker container
# but we want to avoid anyone running multiple instances of bwa index on the
# same reference which will happen if multiple GRIDSS jobs get run in parallel
if [[ ! -f "$REFERENCE.bwt" ]] ; then
	echo "Missing bwa index for $REFERENCE. Attempting to create" 1>&2
	docker run \
		-v "$(dirname $(readlink -f $REFERENCE)):/data/reference/" \
		--rm \
		--entrypoint bwa \
		$CONTAINER \
		index "/data/reference/$(basename $REFERENCE)"
	if [[ ! -f "$REFERENCE.bwt" ]] ; then
		echo "Could not generate $REFERENCE.bwt. Create a bwa index (using \"bwa index $REFERENCE\") or symlink the index files to the expected file names." 1>&2
		exit 1
	fi
fi

# Sanity check coreutils installed on the system
if which readlink 2>&1 >/dev/null ; then
	echo -n
else
	echo "Missing readlink on host. Please install coreutils. Mac users will likely want to run `brew install coreutils`" 1>&2
	exit 1
fi

OUTPUT=${INPUT/.bam/.sv.vcf}
ASSEMBLY=${OUTPUT/.sv.vcf/.gridss.assembly.bam}

RUN_SCRIPT=$OUTDIR/run_gridss_docker_cohort.sh
rm $RUN_SCRIPT 2>/dev/null

grep -vE "^#" $1 | while IFS=',' read -a ALINE ; do
	PATIENT_ID="${ALINE[0]}"
	OUTPUT="${ALINE[1]}"
	ASSEMBLY="${ALINE[2]}"
	if [[ "$OUTPUT" == "" ]] ; then
		OUTPUT=$OUTDIR/$PATIENT_ID.gridss.vcf
	fi
	if [[ "$ASSEMBLY" == "" ]] ; then
		ASSEMBLY=$OUTDIR/$PATIENT_ID.assembly.bam
	fi
	INPUT_ARGS=""
	INPUT_MOUNTS=""
	for ((i = 3; i < ${#ALINE[@]}; i++)) ; do # bams start at the 4th position in the CSV
		BAM="${ALINE[$i]}"
		ORDINAL=$((i - 2)) # shift 
		if [[ ! -f "$BAM" ]] ; then
			echo "Missing input file \"$BAM\"" 1>&2
			exit 1
			echo "DEAD CODE"
		fi
		INPUT_MOUNTS="$INPUT_MOUNTS -v \"$(dirname $(readlink -f $BAM)):/data/bam${ORDINAL}\""
		INPUT_ARGS="$INPUT_ARGS /data/bam${ORDINAL}/$(basename $BAM) "
	done
	PATIENT_SCRIPT=$OUTDIR/run_gridss_docker_patient_${PATIENT_ID}.sh
	cat > $PATIENT_SCRIPT << EOF
#!/bin/sh
docker run \
	--ulimit nofile=$(ulimit -Hn):$(ulimit -Hn) \
	-v "$(dirname $(readlink -f $REFERENCE)):/data/reference/" \
	-v "$(dirname $(readlink -f $BLACKLIST)):/data/blacklist/" \
	-v "$(dirname $(readlink -f $ASSEMBLY)):/data/assembly/" \
	-v "$(dirname $(readlink -f $OUTPUT)):/data/output/" \
	$INPUT_MOUNTS \
	$CONTAINER \
	--workingdir /data/output/ \
	--reference /data/reference/$(basename $REFERENCE) \
	--blacklist /data/blacklist/$(basename $BLACKLIST) \
	--assembly /data/assembly/$(basename $ASSEMBLY) \
	--output /data/output/$(basename $OUTPUT) \
	$INPUT_ARGS \
	2>&1 | tee -a $(readlink -f $OUTDIR)/gridss.$PATIENT_ID.\$HOSTNAME.\$\$.log
EOF
	chmod +x $PATIENT_SCRIPT
	echo "Generated sample script $PATIENT_SCRIPT"
done

#echo "Generated cohort script $RUN_SCRIPT"
