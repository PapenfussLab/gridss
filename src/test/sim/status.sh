#!/bin/bash
# Snippets of usefulness
#
. common.sh

function lockname {
	local LOCK=$(cat $1)
	local NODE=""
	#if [ "$LOCK" != "$(LOCK#unix100)" ] ; then
	NODE=$(qnodes 2>/dev/null | grep -B 4 $LOCK | head -1)
	#fi
	if [ ! -z "$NODE" ] ; then
		NODE="@$NODE"
	fi
	echo $LOCK$NODE
}
TIME_FAIL="Command exited with non-zero status"
# Aligner status
for MD in $DATA_DIR/*.metadata ; do
	cx_load $MD
	ID=$(basename $CX)
	TIME_FILE="${MD%.metadata}.time"
	if [ "$CX_CALLER" != "" ] ; then
		if [[ -f $TIME_FILE ]] && grep "$TIME_FAIL" $TIME_FILE >/dev/null ; then
			echo "$CX_CALLER	$ID	$CX_ALIGNER	${CX_READ_DEPTH}x	${CX_READ_LENGTH}bp	($(tput setaf 6)Execution Failure$(tput sgr0))"
			continue
		elif [ ! -f "$CX_BAM" ] ; then
			echo "$CX_CALLER	$ID	$CX_ALIGNER	${CX_READ_DEPTH}x	${CX_READ_LENGTH}bp	($(tput setaf 45)Missing source bam $CX_BAM$(tput sgr0))"
			continue
		elif [ -f $CX.lock ] ; then
			echo "$CX_CALLER	$ID	$CX_ALIGNER	${CX_READ_DEPTH}x	${CX_READ_LENGTH}bp	($(tput setaf 2)$(lockname $CX.lock)$(tput sgr0))"
			continue
		elif [ -f $CX.vcf -a ! -s $CX.vcf ] ; then
			echo "$CX_CALLER	$ID	$CX_ALIGNER	${CX_READ_DEPTH}x	${CX_READ_LENGTH}bp	($(tput setaf 3)(0 byte VCF)$(tput sgr0))"
			continue
		elif [ ! -f $CX.vcf ] ; then
			echo "$CX_CALLER	$ID	$CX_ALIGNER	${CX_READ_DEPTH}x	${CX_READ_LENGTH}bp	($(tput setaf 1)(Missing VCF)$(tput sgr0))"
			continue
		fi
		# VCF exists happly :)
	elif [ "$CX_ALIGNER" != "" ] ; then
		if [ -f $CX.tmp.bam ] ; then
			if [ -f $CX.lock ] ; then
				echo "$CX_ALIGNER	$ID ($(tput setaf 2)Aligning $(lockname $CX.lock)$(tput sgr0))"
				continue
			else
				echo "$CX_ALIGNER	$ID ($(tput setaf 3)Alignment Failure$(tput sgr0))"
				continue
			fi
		elif [ -f $TIME_FILE ] && grep "$TIME_FAIL" $TIME_FILE >/dev/null ; then
			echo "$CX_ALIGNER	$ID ($(tput setaf 6)Execution Failure$(tput sgr0))"
			continue
		elif [ -f ${TIME_FILE%.time}.sort.time ] && grep "$TIME_FAIL" ${TIME_FILE%.time}.sort.time >/dev/null ; then
			echo "$CX_ALIGNER	$ID ($(tput setaf 6)Sorting execution Failure$(tput sgr0))"
			continue
		elif [ -f $CX.sc.tmp.bam -o -f $CX.sq.tmp.bam ] ; then
			if [ -f $CX.lock ] ; then
				echo "$CX_ALIGNER	$ID ($(tput setaf 2)Sorting $(lockname $CX.lock)$(tput sgr0))"
				continue
			else
				echo "$CX_ALIGNER	$ID ($(tput setaf 3)Sort Failure$(tput sgr0))"
				continue
			fi
		elif [ -f $CX.lock ] ; then
			echo "$CX_ALIGNER	$ID ($(tput setaf 3)$(lockname $CX.lock)$(tput sgr0))"
			continue
		elif [[ ! -f $CX.sam.validation ]] ; then
			echo "$CX_ALIGNER	$ID ($(tput setaf 5)Missing validation file$(tput sgr0))"
			continue
		elif [[ $(grep -c "No errors found" $CX.sam.validation) -eq 0 ]] ; then
			echo "$CX_ALIGNER	$ID ($(tput setaf 5)Validation failure$(tput sgr0))"
			continue
		fi
	elif [ "$CX_READ_SIM" != "" ] ; then
		if [ -f $TIME_FILE ] && grep "$TIME_FAIL" $TIME_FILE >/dev/null ; then
			echo "$CX_READ_SIM	$ID ($(tput setaf 6)Execution Failure$(tput sgr0))"
		elif [ -f $CX.lock ] ; then
			echo "$CX_READ_SIM	$ID ($(tput setaf 2)$(lockname $CX.lock)$(tput sgr0))"
			continue
		elif [ ! -f $CX.2.fq ] ; then
			echo "$CX_READ_SIM	$ID ($(tput setaf 1)(Missing 2.fq)$(tput sgr0))"
			continue
		fi
	elif [ -f $CX.reference.vcf ] ; then
		true # no checks for now
	else
		echo "$(tput setaf 3)Unknown	$ID$(tput sgr0))"
	fi
done
for ID in $(cd $DATA_DIR; ls -1 *.vcf | cut -b 1-32) ; do
	ID=$(basename $ID)
	if [ ! -f $DATA_DIR/$ID.metadata ] ; then
		echo "$ID ($(tput setaf 1)(Missing metadata)$(tput sgr0))"
	fi
done
