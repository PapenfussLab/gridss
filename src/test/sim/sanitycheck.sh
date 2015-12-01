#!/bin/bash
# Snippets of usefulness
#
. common.sh

# Aligner status
for MD in data/*.metadata ; do
	cx_load $MD
	ID=$(basename $CX)
	if [ "$CX_CALLER" != "" ] ; then
		# no generic caller santiy checks yet
		continue
	elif [ "$CX_ALIGNER" != "" ] ; then
		if [ -f $CX.su.bam ] ; then
			# 4 fq lines per read pair = SAM & FQ line counts should match
			FQ_LINES=$(wc -l < $CX_FQ1)
			BAM_LINES=$(samtools flagstat $CX.su.bam | head -1 | cut -d " " -f 1)
			if [ $((FQ_LINES / 2)) -ne "$BAM_LINES" ] ; then
				echo "$(tput setaf 1)$CX_ALIGNER	$ID Expected $((FQ_LINES / 2)) pairs, found $BAM_LINES ($((BAM_LINES - FQ_LINES / 2)))$(tput sgr0)"
			else
				echo "$CX_ALIGNER	$ID $(tput setaf 2)Found $FQ_LINES as expected$(tput sgr0)"
			fi
			if [ -f $CX.sc.bam ] ; then
				if [ ! -f $CX.sq.bam ] ; then
					echo "$(tput setaf 1)$CX_ALIGNER	$ID	Missing $CX.sq.bam$(tput sgr0)"
				else
					if diff <(samtools flagstat $CX.sq.bam) <(samtools flagstat $CX.su.bam) ; then
						echo "$(tput setaf 1)$CX_ALIGNER	$ID	Record mismatch in $CX.sq.bam$(tput sgr0)"
						fi
				fi
				if diff <(samtools flagstat $CX.sc.bam) <(samtools flagstat $CX.su.bam) ; then
					echo "$(tput setaf 1)$CX_ALIGNER	$ID	Record mismatch in $CX.sc.bam$(tput sgr0)"
				fi
			fi
		fi
	elif [ "$CX_READ_SIM" != "" ] ; then
		continue
	elif [ -f $CX.reference.vcf ] ; then
		continue
	else
		continue
	fi
done
