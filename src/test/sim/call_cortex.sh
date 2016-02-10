#!/bin/bash
#
#
. common.sh
CALLER=cortex/1.0.5.14
CORTEX_DIR=$BASE_DIR/tools/cortex/CORTEX_release_v1.0.5.21
# Compilation failed for 1.0.5.21 due to our servers using glibc 2.11
#PATH=/usr/local/bioinfsoftware/cortex/CORTEX_release_v1.0.5.14/bin:$PATH
export PATH=$CORTEX_DIR/bin:$PATH
export PATH=$BASE_DIR/tools/cortex/stampy-1.0.28:$PATH
export PATH=$BASE_DIR/tools/cortex/vcftools_0.1.9/bin:$PATH
export PATH=$CORTEX_DIR/scripts/analyse_variants/needleman_wunsch-0.3.0:$PATH
export PERL5LIB=$CORTEX_DIR/scripts/analyse_variants/bioinf-perl/lib:$PERL5LIB
export PERL5LIB=$CORTEX_DIR/scripts/calling:$PERL5LIB

# Reference precprocessing
# cd $(dirname $CX_REFERENCE)
# stampy.py -G hg19 hg19.fa && stampy.py -g hg19 -H hg19
# cd $CX_REFERENCE.split
# ls -1 $CX_REFERENCE.split/*.fa > $CX_REFERENCE).splitfile_listing_fasta
# cortex_var_31_c1 --kmer_size 31 --mem_height 27 --mem_width 100 --se_list file_listing_fasta --max_read_len 10000 --dump_binary ref.k31.ctx --sample_id REF
# cortex_var_63_c1 --kmer_size 61 --mem_height 27 --mem_width 100 --se_list file_listing_fasta --max_read_len 10000 --dump_binary ref.k61.ctx --sample_id REF

for FQ1 in $DATA_DIR/*.1.generated.fq ; do
	BAM=${FQ1/.1.generated.fq/}
	cx_load $BAM
	# ignore everything that 
	if [[ "$CX_ALIGNER" != "bwamem" ]] ; then
		continue
	fi
	CX_CALLER=$CALLER
	CX_CALLER_ARGS="31,61"
	cx_save
	XC_OUTPUT=$CX.vcf
	# http://cortexassembler.sourceforge.net/cortex_var_user_manual.pdf
	# build ref binaries
	# build stampy hash of genome
	REFSIZE=$(stat -c %s $CX_REFERENCE)
	if [[ -f "$CX_REFERENCE_FA" ]] ; then
		# do we need to adjust for ploidy
		REFSIZE=$(stat -c %s $CX_REFERENCE_FA)
	fi
	XC_SCRIPT="
#rm -rf $CX;
rm -rf $CX/cortex_var/calls $CX/cortex_var/vcfs # don't clobber graph if we have it
mkdir -p $CX 2>/dev/null; cd $CX;
echo \"sample	sample.unpaired.list	sample.paired1.list	sample.paired2.list\" > INDEX
echo $BAM.generated.fq > sample.unpaired.list
echo $BAM.1.generated.fq > sample.paired1.list
echo $BAM.2.generated.fq > sample.paired2.list

perl $CORTEX_DIR/scripts/calling/run_calls.pl \
	--first_kmer 31 \
	--kmer_step 30 \
	--last_kmer 61 \
	--fastaq_index INDEX \
	--auto_cleaning yes \
	--bc yes \
	--pd yes \
	--outdir cortex_var \
	--outvcf cortex.vcf \
	--ploidy 2 \
	--stampy_hash ${CX_REFERENCE/.fa/} \
	--stampy_bin $BASE_DIR/tools/cortex/stampy-1.0.28/stampy.py \
	--list_ref_fasta $CX_REFERENCE.split/file_listing_fasta \
	--refbindir $(dirname $CX_REFERENCE) \
	--genome_size $REFSIZE \
	--qthresh 5 \
	--mem_height 26 --mem_width 100 \
	--vcftools_dir $BASE_DIR/tools/cortex/vcftools_0.1.9/ \
	--do_union yes \
	--ref CoordinatesAndInCalling \
	--workflow independent \
	--logfile log.$(date +%Y%m%d%H%M%S).txt \
	--max_var_len 65537 \
	&& \
	cp $CX/cortex_var/vcfs/cortex.vcf_wk_flow_I_RefCC_FINALcombined_BC_calls_at_all_k.decomp.vcf $CX.vcf && \
	grep -v '#' $CX/cortex_var/vcfs/cortex.vcf_wk_flow_I_RefCC_FINALcombined_PD_calls_at_all_k.decomp.vcf >> $CX.vcf
	
"
	xc_exec
done

