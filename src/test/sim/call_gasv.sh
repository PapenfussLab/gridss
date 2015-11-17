#!/bin/bash
#
# runs GASVPro against bams
#
. common.sh
CALLER=gasv/20140228
JAVA="java -Xmx32g"
for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	if [ $(($CX_READ_FRAGMENT_LENGTH - (2 * $CX_READ_LENGTH + 50 ) )) -le 0 ] ; then
		echo "GASV requires at least 50bp unread within fragments. Skipping Read Length ${CX_READ_LENGTH}bp for ${CX_READ_FRAGMENT_LENGTH}bp fragments."
		continue
	fi
	CHR_PREFIX=""
	if [[ "$(head -c 4 $CX_REFERENCE)" == ">chr" ]] ; then
		CHR_PREFIX="chr"
	fi
	CX_BAM=$BAM
	CX_CALLER=$CALLER
	cx_save
	# use (offset) reference index of chromosomes as custom ordering
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="module add samtools $CALLER ; rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	samtools view -H $BAM | grep -E '^@SQ' | cut -f 2 | cut -b 4- | awk 'BEGIN {I=1} { print \$1 \"\t\" I++ }' > CHROMOSOME_NAMING
	ln -s $BAM $CX/out.bam
	$JAVA -jar \$BAMTOGASV_JAR out.bam -GASVPRO true -LIBRARY_SEPARATED all -OUTPUT_PREFIX out -CHROMOSOME_NAMING CHROMOSOME_NAMING &&\
	$JAVA -jar \$GASV_JAR --output regions --batch out.gasv.in && \
	GASVPro-CC out.gasvpro.in out.gasv.in.clusters && \
	GASVPruneClusters.pl out.gasv.in.clusters.GASVPro.clusters && \
	cp out.gasv.in.clusters.GASVPro.clusters out.gasv.in.clusters.GASVPro.clusters.unconverted && \
	cp out.gasv.in.clusters.GASVPro.clusters.pruned.clusters out.gasv.in.clusters.GASVPro.clusters.pruned.clusters.unconverted && \
	convertClusters out.gasv.in.clusters.GASVPro.clusters.pruned.clusters && \
	convertClusters out.gasv.in.clusters.GASVPro.clusters && \
	$BASE_DIR/gasv2vcf.py < out.gasv.in.clusters.GASVPro.clusters CHROMOSOME_NAMING > $CX.vcf"
	xc_exec
done

