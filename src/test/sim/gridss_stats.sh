#!/bin/bash
#
#
. common.sh

echo "aligner,softclip,readlength,readdepth,fraglength"
for F in $(grep -l gridss $DATA_DIR/*.metadata | cut -b 1-32) ; do
	CX=$F
	VCF=$F.vcf
	if [[ -f $F.vcf ]] ; then
		cx_load $F
		
	fi
done
