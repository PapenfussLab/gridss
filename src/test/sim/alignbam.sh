#!/bin/bash
#
# Aligns paired fastq reads
#
. common.sh

PATH=$PATH:$BASE_DIR/tools/variationhunter/mrfast-2.6.0.11

function align_subread {
	cx_load $1
	CX_FQ1=$1
	CX_FQ2=$2
	CX_ALIGNER=subread
	CX_ALIGNER_SOFTCLIP=1
	cx_save
	INDEX=$CX_REFERENCE.subread
	if [ ! -f $INDEX.00.b.array ] ; then
		echo "Unable to find subread index for $CX_REFERENCE"
		echo subread-buildindex -o $CX_REFERENCE.subread $CX_REFERENCE 
		exit 1
	fi
	if [[ "$1" == *.gz ]] ; then
		echo "subread does not handle .gz compressed input, skipping $1"
		return
	fi 
	XC_MULTICORE=1
	XC_OUTPUT=$CX.su.bam
	XC_SCRIPT="module add subread
	subread-align -T $XC_CORES -i $INDEX -r $1 -R $2 | \
	AddOrReplaceReadGroups I=/dev/stdin O=$CX.tmp.bam $READ_GROUP_PARAMS && \
	mv $CX.tmp.bam $XC_OUTPUT
	"
	xc_exec
}
function align_bwa {
	cx_load $1
	CX_FQ1=$1
	CX_FQ2=$2
	CX_ALIGNER=bwamem
	CX_ALIGNER_FLAGS=$3
	CX_ALIGNER_SOFTCLIP=1
	cx_save
	INDEX=$(ls -1 ${CX_REFERENCE/.fa/}*.bwt)
	INDEX=${INDEX/.bwt/}
	if [ ! -f $INDEX.bwt ] ; then
		echo "Unable to find bwa index $INDEX for $CX_REFERENCE"
		echo bwa index $CX_REFERENCE
		exit 1
	fi
	XC_MULTICORE=1
	XC_OUTPUT=$CX.su.bam
	XC_SCRIPT="
	bwa mem -t $XC_CORES $3 $INDEX $1 $2 \
	| AddOrReplaceReadGroups I=/dev/stdin O=$CX.tmp.bam $READ_GROUP_PARAMS && \
	mv $CX.tmp.bam $XC_OUTPUT
	"
	xc_exec
}
function align_bowtie2 {
	cx_load $1
	CX_FQ1=$1
	CX_FQ2=$2
	CX_ALIGNER=bowtie2
	if [ "$3" == "noSoftClip" ] ; then
		CX_ALIGNER_FLAGS="--end-to-end"
		CX_ALIGNER_SOFTCLIP=0
	else
		CX_ALIGNER_FLAGS="--local"
		CX_ALIGNER_SOFTCLIP=1
		
	fi
	cx_save
	#INDEX=$(ls -1 ${CX_REFERENCE/.fa/}*.rev.1.bt2)
	#INDEX=${INDEX/.rev.1.bt2/}
	INDEX=$CX_REFERENCE
	if [[ ! -f $INDEX.rev.1.bt2 ]] ; then
		echo "Unable to find bowtie2 index for $CX_REFERENCE"
		echo bowtie2-build $CX_REFERENCE $CX_REFERENCE
		exit 1
	fi
	XC_MULTICORE=1
	XC_OUTPUT=$CX.su.bam
	XC_SCRIPT="
	bowtie2 --threads $XC_CORES --mm $CX_ALIGNER_FLAGS -x $INDEX -1 $1 -2 $2 \
	| AddOrReplaceReadGroups I=/dev/stdin O=$CX.tmp.bam $READ_GROUP_PARAMS && \
	mv $CX.tmp.bam $XC_OUTPUT
	"
	xc_exec
}
function align_novoalign {
	cx_load $1
	CX_FQ1=$1
	CX_FQ2=$2
	CX_ALIGNER=novoalign
	if [ "$3" == "noSoftClip" ] ; then
		CX_ALIGNER_FLAGS="-o FullNW"
		CX_ALIGNER_SOFTCLIP=0
	else
		CX_ALIGNER_FLAGS="-o Softclip"
		CX_ALIGNER_SOFTCLIP=1
	fi
	cx_save
	INDEX=$(ls -1 ${CX_REFERENCE/.fa/}*.novoindex)
	if [ ! -f "$INDEX" ] ; then
		echo "Unable to find novoalign index for $CX_REFERENCE"
		novoindex $CX_REFERENCE.novoindex $CX_REFERENCE 
		exit 1
	fi
	if [ $(( $(stat -c%s $INDEX) / 1000000 )) -gt $(( $MAX_MEMORY  - 1024 )) ] ; then
		echo "$(basename $CX): Not scheduling novoalign as index $INDEX requires more memory than is available: $(($(stat -c%s $INDEX) / 1000000))Mb required, $(( $MAX_MEMORY  - 1024 ))Mb free"
		return
	fi
	XC_MULTICORE=1
	XC_OUTPUT=$CX.su.bam
	XC_SCRIPT="
	novoalign -c $XC_CORES -o SAM $CX_ALIGNER_FLAGS -F STDFQ -i PE $CX_READ_FRAGMENT_LENGTH,$CX_READ_FRAGMENT_STDDEV \
		-d $INDEX -f $1 $2 \
	| AddOrReplaceReadGroups I=/dev/stdin O=$CX.tmp.bam $READ_GROUP_PARAMS && \
	mv $CX.tmp.bam $XC_OUTPUT
	"
	xc_exec
}

function align_mrfast {
	cx_load $1
	CX_ALIGNER_MIN_CONCORD=$((CX_READ_FRAGMENT_LENGTH - 4 * CX_READ_FRAGMENT_STDDEV ))
	CX_ALIGNER_MAX_CONCORD=$((CX_READ_FRAGMENT_LENGTH + 4 * CX_READ_FRAGMENT_STDDEV ))
	CX_FQ1=$1
	CX_FQ2=$2
	CX_ALIGNER=mrfast
	# TODO: confirm mrfast does not softclip
	CX_ALIGNER_SOFTCLIP=0
	cx_save
	INDEX=$(ls -1 $CX_REFERENCE.index)
	if [ ! -f "$INDEX" ] ; then
		echo "Fatal: Unable to find mrfast index for $CX_REFERENCE"
		echo mrfast --index $CX_REFERENCE
		exit 1
	fi
	XC_OUTPUT=$CX.su.bam
	# 3 stddev = concordant
	MIN=$((CX_READ_FRAGMENT_LENGTH - 3 * CX_READ_FRAGMENT_STDDEV ))
	MAX=$((CX_READ_FRAGMENT_LENGTH + 3 * CX_READ_FRAGMENT_STDDEV ))
	XC_SCRIPT="
	mrfast \
		--search $CX_REFERENCE \
		--pe \
		--discordant-vh \
		--seq1 $CX_FQ1 \
		--seq2 $CX_FQ2 \
		-o $CX.output \
		-u $CX.unmapped \
		--min $CX_ALIGNER_MIN_CONCORD \
		--max $CX_ALIGNER_MAX_CONCORD && \
	AddOrReplaceReadGroups I=<(samtools view -ShT $CX_REFERENCE $CX.output_BEST.sam) O=$CX.tmp.bam $READ_GROUP_PARAMS && \
	mv $CX.tmp.bam $XC_OUTPUT
	"
	xc_exec
}
# Align reads
for FQ1 in $(ls -1 $DATA_DIR/*.1.fastq.gz $DATA_DIR/*.1.fq 2>/dev/null) ; do
	FQ2=${FQ1/.1./.2.}
	cx_load $FQ1
	if [ ! -z "$CX_DOWNSAMPLE_FROM" ] ; then
		# skip downsampled bams as they can be efficiently generated
		# by downsamplebam.sh
		continue
	fi
	if [ -f $CX.lock ] ; then
		echo "$FQ1 still being generated: skipping"
		continue
	fi
	#FQ1="<( bamToFastq -i $CX_BAM -fq /dev/stdout -fq2 /dev/null )"
	#FQ2="<( bamToFastq -i $CX_BAM -fq /dev/null -fq2 /dev/stdout )"
	READ_GROUP_PARAMS="RGPL=illumina RGPU=sim RGSM=$(basename $CX_REFERENCE_VCF) RGLB=$(basename ${CX_REFERENCE_VCF})_rl${CX_READ_LENGTH}_f${CX_READ_FRAGMENT_LENGTH} RGCN=sim${CX_READ_SIM}"
	if [ "$1" != "" ] ; then
		align_$1 $FQ1 $FQ2
		if [ "$1" == "novoalign" -o "$1" == "bowtie2" ] ; then
			align_$1 $FQ1 $FQ2 "noSoftClip"
		fi
	else
		align_bwa $FQ1 $FQ2
		if [[ $FULL_MATRIX -eq 1 ]] ; then
			align_bwa $FQ1 $FQ2 "-M -Y"
			align_novoalign $FQ1 $FQ2
			#align_subread $FQ1 $FQ2 # exclude subread until MAPQ & SC issues resolved
			align_bowtie2 $FQ1 $FQ2
			
			#align_novoalign $FQ1 $FQ2 "noSoftClip"
			#align_bowtie2 $FQ1 $FQ2 "noSoftClip"
			#align_mrfast $FQ1 $FQ2
		fi
	fi
done

