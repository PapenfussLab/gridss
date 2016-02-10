#!/bin/bash
#
# runs breakdancer against bams
#
. common.sh
CALLER=meerkat/0.185
MEERKAT_DIR=$BASE_DIR/tools/meerkat
export LD_LIBRARY_PATH=$MEERKAT_DIR/src/mybamtools/lib:$LD_LIBRARY_PATH
export PATH=$MEERKAT_DIR/bin:$PATH
export PATH=/usr/local/bioinfsoftware/samtools/samtools-0.1.16:$PATH
export PATH=/usr/local/bioinfsoftware/bwa/bwa-0.6.2/bin:$PATH
export PATH=/usr/local/bioinfsoftware/blast/blast-2.2.18/bin:$PATH

for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	if [[ ! -d $CX_REFERENCE.meerkat ]] ; then
		# create isolated reference directory as required my meerkat
		mkdir $CX_REFERENCE.meerkat
		ln -s $CX_REFERENCE $CX_REFERENCE.meerkat/$(basename $CX_REFERENCE)
	fi
	CX_BAM=$BAM
	CX_CALLER=$CALLER
	CX_CALLER_FLAGS=""
	cx_save
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	mkdir input
	ln -s $CX_BAM input/input.bam
	ln -s $CX_BAM.bai input/input.bam.bai
	perl $MEERKAT_DIR/scripts/pre_process.pl \
		-b input/input.bam \
		-I $CX_REFERENCE \
		-A $CX_REFERENCE.fai \
		
	perl $MEERKAT_DIR/scripts/meerkat.pl \
		-b input/input.bam \
		-u 1 \
		-Q 10 \
		-F $CX_REFERENCE.meerkat \
		
	perl $MEERKAT_DIR/scripts/mechanism.pl \
		-b input/input.bam \
		-R $CX_REFERENCE.ucsc.table.rmsk.txt \
		
	perl $MEERKAT_DIR/scripts/meerkat2vcf.pl \
		-i input/input.variants \
		-o out.vcf \
		-H $BASE_DIR/meerkat.vcf \
		-F $CX_REFERENCE.meerkat \
		
	cp out.vcf $XC_OUTPUT
		
	"
	xc_exec
done

