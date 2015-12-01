#!/bin/bash
. common.sh

echo "caller	readDepth	readLength	fragmentLength	aligner	softclipped	Huser	user	Hreal	real	Hsys	sys	HmaxResMem	maxResMem	HavgTotMem	avgTotMem	HI	I	HO	O	Hexit	exit"
for MD in $DATA_DIR/*.metadata ; do
	cx_load $MD
	ID=$(basename $CX)
	if [ "$CX_CALLER" != "" ] ; then
		if [ -f $CX.lock ] ; then
			echo "$CX_CALLER	$ID ($(tput setaf 2)$(cat $CX.lock)$(tput sgr0))" 1>&2
			continue
		elif [ -f $CX.vcf -a ! -s $CX.vcf ] ; then
			echo "$CX_CALLER	$ID ($(tput setaf 3)(0 byte VCF)$(tput sgr0))" 1>&2
			continue
		elif [ ! -f $CX.vcf ] ; then
			echo "$CX_CALLER	$ID ($(tput setaf 1)(Missing VCF)$(tput sgr0))" 1>&2
			continue
		elif [ ! -f $CX.time ] ; then
			# check for POSIX time info in stderr
			if [[ ! "$(tail -3 $CX.stderr | grep real)" =~ ([0-9]+)m([0-9]+)\.[0-9]*s ]] ; then
				echo "$CX_CALLER	$ID ($(tput setaf 1)(Missing .time and no timing information in stderr)$(tput sgr0))" 1>&2
				continue
			fi
			REAL_TIME=$(( ${BASH_REMATCH[1]} * 60 + ${BASH_REMATCH[2]} ))
			[[ ! "$(tail -3 $CX.stderr | grep user)" =~ ([0-9]+)m([0-9]+)\.[0-9]*s ]]
			USER_TIME=$(( ${BASH_REMATCH[1]} * 60 + ${BASH_REMATCH[2]} ))
			[[ ! "$(tail -3 $CX.stderr | grep sys)" =~ ([0-9]+)m([0-9]+)\.[0-9]*s ]]
			SYS_TIME=$(( ${BASH_REMATCH[1]} * 60 + ${BASH_REMATCH[2]} ))
			RUNTIME="user	$USER_TIME	elapsed	$REAL_TIME	sys	$SYS_TIME	maxResMemKb	0	avgTotMemKb	0	I	0	O	0	exit	0"
		elif grep "non-zero" $CX.time >/dev/null ; then
			echo "$CX_CALLER	$ID ($(tput setaf 1)(Execution failure)$(tput sgr0))" 1>&2
			continue
		else
			RUNTIME=$(cat $CX.time)
		fi
		echo "$CX_CALLER	$CX_READ_DEPTH	$CX_READ_LENGTH	$CX_READ_FRAGMENT_LENGTH	$CX_ALIGNER	$CX_ALIGNER_SOFTCLIP	$RUNTIME"
	fi
done
