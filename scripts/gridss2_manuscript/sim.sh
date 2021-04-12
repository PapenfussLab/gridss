#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=8g
#SBATCH --time=2:00:00
#SBATCH -p regular
#SBATCH --output=bench.%J.out.%x.log
#SBATCH --error=bench.%J.err.%x.log
#SBATCH --array=1-3004
source ~/miniconda3/etc/profile.d/conda.sh
# conda create -n gridss2_benchmark_sim manta=1.6.0 gridss=2.11.1 svaba=1.1.0 bwa=0.7.17 samtools=1.11 bioconductor-structuralvariantannotation art=2016.06.05
# novobreak conda doesn't work :( novobreak=1.1.3rc
conda activate gridss2_benchmark_sim
threads=4
depth=60 # manta fails at depth=30 due to insufficient RPs to infer library fragment size distribution
flank_size=10000
art_args="--noALN --paired --seqSys HSXt -l 150 -m 500 -s 100 -rs 1"
full_ref=~/dev/ref/hg19.fa
base_dir=$(readlink -f publicdata/sim/)
gen_dir=$base_dir/gen
mkdir -p $gen_dir
if [[ -f ~/dev/gridss/target/github_package/gridss-2.11.1-gridss-jar-with-dependencies.jar ]] ; then
	export GRIDSS_JAR=~/dev/gridss/target/github_package/gridss-2.11.1-gridss-jar-with-dependencies.jar
fi
if [[ ! -f "$GRIDSS_JAR" ]] ; then
	echo "Missing GRIDSS_JAR"
	exit 1
fi
t_all_sim_fasta=$base_dir/all_runs.t.fa
n_all_sim_fasta=$base_dir/all_runs.n.fa
if [[ ! -f $n_all_sim_fasta ]] ; then
	echo "Missing $n_all_sim_fasta" 1>&2
	echo "Generate with gen_sim.R" 1>&2
	exit 1
fi
if [[ ! -f $t_all_sim_fasta.fai ]] ; then 
	samtools faidx $t_all_sim_fasta
fi
if [[ ! -f $n_all_sim_fasta.fai ]] ; then 
	samtools faidx $n_all_sim_fasta --rcount
fi
if [[ ! -d $base_dir/nb_distribution ]] ; then
	git clone https://github.com/czc/nb_distribution.git $base_dir/nb_distribution
fi
#novobreak is really slow against a full reference: we'll just use chr8
ref=$base_dir/$(basename $full_ref).chr8.fa
if [[ ! -f $ref.gridsscache ]] ; then
	samtools faidx $full_ref chr8 > $ref
	samtools faidx $ref
	gridss.sh -j $GRIDSS_JAR -s setupreference -r $ref -o ignored.vcf $ref
fi
ref_reads_bam=$gen_dir/refreads.bam
ref_reads_fasta=$gen_dir/refreads.fa
ref_reads_fastq=$gen_dir/refreads.fq
if [[ ! -f $ref_reads_bam ]] ; then
	# hg19 has no gap (N bases) here
	samtools faidx $ref chr8:1000000-1999999 > $ref_reads_fasta
	# SVABA requires 2M reads
	art_illumina $art_args -c 1000000 --id ref_ --in $ref_reads_fasta -o $ref_reads_fastq.
	bwa mem -t $(nproc) $ref $ref_reads_fastq.1.fq $ref_reads_fastq.2.fq | samtools sort -@ $(nproc) -O BAM -o $ref_reads_bam -
	samtools index $ref_reads_bam
	gzip $gen_dir/refreads.fq.*
fi
# $1 run index
process_run() {
	n=$1
	run=run${n}_
	run_dir=$gen_dir/$run
	tsrc_fa=$run_dir/t.${run}.fa
	nsrc_fa=$run_dir/n.${run}.fa
	tfq1=$run_dir/t.$run.${depth}x.1.fq
	tfq2=$run_dir/t.$run.${depth}x.2.fq
	nfq1=$run_dir/n.$run.${depth}x.1.fq
	nfq2=$run_dir/n.$run.${depth}x.2.fq
	tbam_raw=$run_dir/t.$run.targeted.bam
	nbam_raw=$run_dir/n.$run.targeted.bam
	tbam_merged=$run_dir/t.$run.merged.bam
	nbam_merged=$run_dir/n.$run.merged.bam
	mkdir -p $run_dir
	cd $run_dir
	if [[ ! -f $nfq2 ]] ; then
		samtools faidx $t_all_sim_fasta run_${n}_ > $tsrc_fa
		samtools faidx $n_all_sim_fasta anchor_run_${n}_ centro_run_${n}_  > $nsrc_fa
		art_illumina $art_args --fcov $depth --id t${n}_ --in $tsrc_fa -o $run_dir/t.$run.${depth}x.
		art_illumina $art_args --fcov $depth --id n${n}_ --in $nsrc_fa -o $run_dir/n.$run.${depth}x.
	fi
	if [[ ! -f $nbam_merged ]] ; then
		bwa mem -t $threads $ref $tfq1 $tfq2 | samtools sort -@ $threads -O BAM -o $tbam_raw -
		samtools index $tbam_raw
		samtools merge $tbam_merged $ref_reads_bam $tbam_raw
		samtools index $tbam_merged
		bwa mem -t $threads $ref $nfq1 $nfq2 | samtools sort -@ $threads -O BAM -o $nbam_raw -
		samtools index $nbam_raw
		samtools merge $nbam_merged $ref_reads_bam $nbam_raw
		samtools index $nbam_merged
	fi
	for caller in gridss manta svaba novobreak ; do
		echo "Processing $caller $n"
		caller_dir=$run_dir/$caller
		mkdir -p $caller_dir
		vcf=$run_dir/$run.$caller.vcf
		if [[ ! -f $vcf ]] ; then
			cd $caller_dir
			if [[ $caller == "gridss" ]] ; then
				echo "chunkSize=1000000000" > gridss.config
				gridss.sh -c gridss.config -t $threads --jvmheap 2g -o $(basename $vcf) -r $ref -j $GRIDSS_JAR ../$(basename $nbam_merged) ../$(basename $tbam_merged)
				gridss_somatic_filter.R -i $(basename $vcf) -o $(basename $vcf).somatic.vcf --fulloutput $(basename $vcf).full.vcf -s ~/dev/gridss/scripts -c ~/dev/gridss/scripts -t 2 -n 1
				gunzip -c $(basename $vcf).somatic.vcf.bgz > $vcf
			elif [[ $caller == "manta" ]] ; then
				configManta.py --normalBam=../$(basename $nbam_merged) --tumourBam=../$(basename $tbam_merged) --runDir=. --referenceFasta=$ref
				./runWorkflow.py -j $threads
				cp $caller_dir/workspace/svHyGen/somaticSV.0000.vcf $vcf
			elif [[ $caller == "svaba" ]] ; then
				svaba run -t ../$(basename $tbam_merged) -n ../$(basename $nbam_merged) -p $threads -a $run -G $ref
				cp $caller_dir/run1_.svaba.somatic.sv.vcf $vcf
			elif [[ $caller == "novobreak" ]] ; then
				oldpath=$PATH
				ln -s ../$(basename $tbam_raw)
				ln -s ../$(basename $nbam_raw)
				export PATH=$base_dir/nb_distribution/:$PATH
				run_novoBreak.sh $base_dir/nb_distribution $ref $(basename $tbam_raw) $(basename $nbam_raw) $threads
				export PATH=$oldpath
				cp $caller_dir/novoBreak.pass.flt.vcf $vcf
			else
				echo "ERROR"
				exit 1
			fi
			cd $run_dir
		fi
	done
}
if [[ "$1" != "" ]] ; then
	process_run $1
elif [[ "$SLURM_ARRAY_TASK_ID" != "" ]] ; then
	process_run $1
else
	for n in {1..3004..1024} ; do
		process_run $n
	done
fi

