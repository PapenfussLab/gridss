#!/bin/bash
#
# Generates simulated read data sets based on read error profiles
#
. common.sh

# ART
# pIRS http://www.ncbi.nlm.nih.gov/pubmed/22508794

ART_DIR=$BASE_DIR/tools/art
PATH=$PATH:$ART_DIR/art_illumina_dir/src:$ART_DIR/Linux64

# Generate simulated reads for each reference VCF
for VCF in $(ls_reference_vcf) ; do
	echo "Simulating reads using ART default profiles"
	for READ_LENGTH in $READ_LENGTHS; do
		case $READ_LENGTH in
		36 | 44 | 50 | 75)
			QUAL_FILE_1=$ART_DIR/art_illumina_dir/Illumina_GAII_profiles/EmpR${READ_LENGTH}R1.txt
			QUAL_FILE_2=$ART_DIR/art_illumina_dir/Illumina_GAII_profiles/EmpR${READ_LENGTH}R2.txt
			;;
		100 ) 
			QUAL_FILE_1=$ART_DIR/art_illumina_dir/Illumina_GAII_profiles/Emp${READ_LENGTH}R1.txt
			QUAL_FILE_2=$ART_DIR/art_illumina_dir/Illumina_GAII_profiles/Emp${READ_LENGTH}R2.txt
			;;
		150 | 200 | 250 | * )
			QUAL_FILE_1=$ART_DIR/art_illumina_dir/Illumina_GAII_profiles/EmpMiSeq250R1.txt
			QUAL_FILE_2=$ART_DIR/art_illumina_dir/Illumina_GAII_profiles/EmpMiSeq250R2.txt
			;;
		esac
# temp hack: are pindel & crest bad due to quality filtering?
QUAL_FILE_1=$ART_DIR/art_illumina_dir/Illumina_GAII_profiles/EmpMiSeq250R1.txt
QUAL_FILE_2=$ART_DIR/art_illumina_dir/Illumina_GAII_profiles/EmpMiSeq250R2.txt
		if [ ! -f $QUAL_FILE_1 ] ; then
			echo "Missing art base quality file $QUAL_FILE_1" 1>&2
			continue
		fi
		if [ ! -f $QUAL_FILE_1 ] ; then
			echo "Missing art base quality file $QUAL_FILE_1" 1>&2
			continue
		fi
		for FRAG_SIZE in $FRAGMENT_SIZE ; do
			# only do other read lengths for 300bp fragment size
			if [[ $FRAG_SIZE -eq 300 || $READ_LENGTH -eq 100 || $FULL_MATRIX -eq 1 ]] ; then
				cx_load $VCF
				CX_REFERENCE_VCF=$VCF
				CX_REFERENCE_FA=${VCF%vcf}fa
				CX_READ_LENGTH=$READ_LENGTH
				CX_READ_FRAGMENT_LENGTH=$FRAG_SIZE
				CX_READ_FRAGMENT_STDDEV=$(( CX_READ_FRAGMENT_LENGTH / 10))
				CX_READ_DEPTH=$STARTING_DEPTH
				CX_READ_SIM="art"
				CX_READ_SIM_QUAL_FILE=$QUAL_FILE_1
				#CX_READ_SIM_SEED=12345
				cx_save
				XC_OUTPUT=$CX.1.fq
				XC_SCRIPT="
				rm $CX.1.fq $CX.2.fq 2>/dev/null
				art_illumina \
					--paired \
					--in $CX_REFERENCE_FA \
					--out $CX. \
					--noALN \
					--len $READ_LENGTH \
					--mflen $FRAG_SIZE \
					--sdev $CX_READ_FRAGMENT_STDDEV \
					--fcov $((CX_READ_DEPTH / 2)) \
					--rndSeed 0 \
					--qprof1 $QUAL_FILE_1 \
					--qprof2 $QUAL_FILE_2 \
					--id art
				fastqc --threads 2 --nogroup -o $DATA_DIR -f fastq $CX.1.fq $CX.2.fq
				"
				xc_exec
			fi
		done
	done
done
