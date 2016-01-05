#!/bin/bash
#
# runs breakdancer against bams
#
unset PERL5LIB # temp hack until environment gets fixed
. common.sh
CALLER=crest/0.0.1
export PATH=/usr/local/bioinfsoftware/blat/blat_34x12/bin:$PATH 
BLAT_PORT=24513
BLAT_POLL_INTERVAL=10
for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=crest
	CX_CALLER_ARGS=
	cx_save
	if [ ! -f $CX_REFERENCE.2bit ] ; then
		echo "Missing $CX_REFERENCE.2bit - unable run CREST"
		exit 1
	fi
	BLAT_TMP_FILE=/tmp/crest.blat.$BLAT_PORT.$(basename $CX).lock
	XC_OUTPUT=$CX.vcf
	XC_TRAP="rm $BLAT_TMP_FILE ; if ls /tmp/crest.blat.$BLAT_PORT.*.lock 2>/dev/null ; then echo 'Additional instance using blat server' ; else gfServer stop localhost $BLAT_PORT 2>/dev/null; pkill -s 0 gfServer; echo 'Blat server terminated' ; fi "
	XC_SCRIPT="module add cap3 ucsc-tools $CALLER ; rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	ln -s $BAM $CX/in.bam
	ln -s $BAM.bai $CX/in.bam.bai
	echo > $BLAT_TMP_FILE
	gfServer start localhost $BLAT_PORT $CX_REFERENCE.2bit & 
	while
		echo Waiting for ${BLAT_POLL_INTERVAL}s for BLAT server to accept requests 1>&2
		sleep $BLAT_POLL_INTERVAL
		gfClient localhost $BLAT_PORT $CX_REFERENCE.2bit <( echo -e '>\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' ) /dev/null
		[[ \$? != 0 ]]
	do
		:
	done
	echo Success connecting to BLAT server 1>&2
	# output file is in working directory so we need to move to the data directory before running CREST
	extractSClip.pl -i in.bam --ref_genome $CX_REFERENCE
	CREST.pl $CX_CALLER_ARGS -l $CX_READ_LENGTH -f in.bam.cover -d in.bam --ref_genome $CX_REFERENCE -t $CX_REFERENCE.2bit --blatserver localhost --blatport $BLAT_PORT && \
	$BASE_DIR/crest2vcf.py < in.bam.predSV.txt > $CX.vcf
	"
	xc_exec
done

