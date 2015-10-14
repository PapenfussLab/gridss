#!/bin/bash
#
# Runs SVDetect against BAMs
#
. common.sh

SVDETECT_DIR=$BASE_DIR/tools/svdetect/SVDetect_r0.8b
PERL5LIB=$PERL5LIB:$SVDETECT_DIR/lib
PATH=$PATH:$SVDETECT_DIR/bin

for BAM in $DATA_DIR/*.sq.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=svdetect
	cx_save
	cat > $CX.conf << EOF
<general>
input_format=sam/bam
sv_type=all
mates_orientation=RF
read1_length=$CX_READ_LENGTH
read2_length=$CX_READ_LENGTH
mates_file=$CX.discordant.sam
cmap_file=$CX.len
output_dir=$CX
tmp_dir=$CX/tmp
num_threads=8
</general>
<detection>
split_mate_file=1
window_size=5000
step_length=1000
</detection>
<filtering>
split_link_file=1
nb_pairs_threshold=3
strand_filtering=1
</filtering>
EOF
	# generate the chromosome length file from the fa index file
	# as it already has the required info
 	awk ' { print NR "\t" $1 "\t" $2 }' < $CX_REFERENCE.fai > $CX.len
	XC_OUTPUT=$CX.discordant.sam.links.filtered.sv.txt
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	mkdir $CX/tmp
	# remove proper pairs
	# remove unmapped
	# remove mate unmapped
	# remove secondary alignments
	# remove supplementary alignments
	# require paired read
	samtools view -hF 2 $CX_BAM \
	| samtools view -ShF 4 - \
	| samtools view -ShF 8 - \
	| samtools view -ShF 256 - \
	| samtools view -ShF 2048 - \
	| samtools view -Shf 1 - > $CX.discordant.sam
	SVDetect linking -conf $CX.conf
	SVDetect filtering -conf $CX.conf
	SVDetect links2SV -conf $CX.conf
	"
	# SVDetect 0.8b doesn't seem to run:
	# Use of uninitialized value $strand1 in string eq at /usr/local/bioinf/bin/SVDetect line 499, <MATES> line 1.
	# Use of uninitialized value $strand2 in string eq at /usr/local/bioinf/bin/SVDetect line 500, <MATES> line 1.
	# Use of uninitialized value $chr_read1 in exists at /usr/local/bioinf/bin/SVDetect line 542, <MATES> line 1.
	xc_exec
done

