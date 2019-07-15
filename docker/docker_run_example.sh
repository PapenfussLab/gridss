#!/bin/bash

INPUT=../example/chr12.1527326.DEL1024.bam
REFERENCE=../hg19.fa
BLACKLIST=../example/wgEncodeDacMapabilityConsensusExcludable.bed
OUTDIR=gridss.example.docker
OUTPUT=$OUTDIR/$(basename ${INPUT/.bam/}.gridss.vcf.gz)
mkdir -p $OUTDIR
CMD=$(echo docker run \
	--ulimit nofile=$(ulimit -Hn):$(ulimit -Hn) \
	-v "$(dirname $(readlink -f $REFERENCE)):/data/reference/" \
	-v "$(dirname $(readlink -f $BLACKLIST)):/data/blacklist/" \
	-v "$(dirname $(readlink -f $INPUT)):/data/input/" \
	-v "$(dirname $(readlink -f $OUTPUT)):/data/output/" \
	gridss/gridss:latest \
	-r /data/reference/$(basename $REFERENCE) \
	-b /data/blacklist/$(basename $BLACKLIST) \
	-w /data/output/ \
	-o /data/output/$(basename $OUTPUT) \
	--assembly /data/output/${OUTPUT/.vcf.gz/}.assembly.bam \
	/data/input/$(basename $INPUT)
)
echo Linux:
echo $CMD
echo
echo WSL:
echo ${CMD//\/mnt\/d/d:/}
echo
