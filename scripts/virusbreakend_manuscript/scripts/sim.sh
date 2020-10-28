#!/bin/bash
#conda create -n vifi bwa=0.7.17 python=2.7 pysam=0.15.2 samtools=1.9 hmmer
export PATH=~/dev/gridss/:$PATH
export GRIDSS_JAR=~/projects/virusbreakend/gridss-2.10.1-gridss-jar-with-dependencies.jar
hpv=MT783410.1
hbv=LC500247.1
chr1_len=248387497
n=0
pos=0
flank=50000
breakoffset=10
virus_lengh=1000
readlen=150
virusbreakenddb=~/projects/virusbreakend/virusbreakenddb
ref=~/projects/reference_genomes/human/hg19.fa
if [[ ! -f gen/chm13.bam.gridss.working/chm13.bam.insert_size_metrics ]] ; then
	# GRIDSS uses WGS metrics which are going to be wrong for the small simulation - create a baseline
	wgsim -s 1 -e 0.005 -r 0 -R 0 -X 0 -N 20000000 -1 $readlen -2 $readlen chm13.draft_v1.0.fasta gen/chm13.1.fq gen/chm13.2.fq
	bwa mem -t $(nproc) $ref gen/chm13.1.fq gen/chm13.2.fq | samtools sort -@ $(nproc) -o gen/chm13.bam -O BAM -
	gridss.sh -r $ref -s preprocess gen/chm13.bam -w gen/
fi
while [[ $n -lt 248 ]] ; do
	n=$(($n + 1))
	pos=$(($n * 1000000))
	virus_pos=$(($n * 8))
	file=$PWD/gen/chr1_${pos}.fa
	if [[ ! -f $file ]] ; then
		echo ">$file" >> $file
		samtools faidx chm13.draft_v1.0.fasta chr1:$(($pos - $flank))-$(($pos)) | grep -v ">" | tr -d '\n' >> $file
		samtools faidx $hbv.fa $hbv:$virus_pos-$(($virus_pos + $virus_lengh)) | grep -v ">" | tr -d '\n' >> $file
		samtools faidx chm13.draft_v1.0.fasta chr1:$(($pos + $breakoffset))-$(($pos + $flank + $breakoffset)) | grep -v ">" | tr -d '\n' >> $file
		echo "Generated $n"
	fi
	reflen=$(( 2 * $flank + $virus_lengh ))
	pairbases=$(( 2 * $readlen ))
	for depth in 5 10 15 30 60 ; do
		if [[ ! -f $file.${depth}x.1.fq ]] ; then
			wgsim -S 1 -r 0 -R 0 -X 0 -N $(( $depth * $reflen / $pairbases )) -1 150 -2 150 $file $file.${depth}x.1.fq $file.${depth}x.2.fq
			cat ecoli_1.fq >> $file.${depth}x.1.fq
			cat ecoli_2.fq >> $file.${depth}x.2.fq
			sed -i 's/222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222/FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF2/' $file.${depth}x.1.fq
			sed -i 's/222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222/FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF2/' $file.${depth}x.2.fq
		fi
		vcf=gen/virusbreakend_${n}_${depth}x.vcf
		bam=$file.${depth}x.bam
		if [[ ! -f $bam ]]; then
				cat > gen/slurm_bwa_${n}_${depth}x.sh << EOF
#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=10g
#SBATCH --time=12:00:00
#SBATCH -p regular
#SBATCH --output=log.bwa.%x.out
#SBATCH --error=log.bwa.%x.err
. ~/conda_crest/etc/profile.d/conda.sh
conda activate virusbreakend
cd $PWD
cd gen
if [[ ! -f $bam ]] ; then
	echo "Creating $bam"
	bwa mem -t 4 $ref $file.${depth}x.1.fq $file.${depth}x.2.fq | samtools sort -@ 4 -o $bam -O BAM -
	samtools index $bam
fi
EOF
		fi
		if [[ ! -f $vcf.virusbreakend.working/$(basename $vcf).kraken2.report.viral.extracted.txt ]] ; then
			metrics_prefix=$vcf.virusbreakend.working/adjusted/$(basename $bam).viral.bam.gridss.working/$(basename $bam).viral.bam
			mkdir -p $(dirname $metrics_prefix)
			for metric_suffix in cigar_metrics idsv_metrics insert_size_metrics mapq_metrics tag_metrics ; do
				if [[ ! -f $metrics_prefix.$metric_suffix ]] ; then
					cp gen/chm13.bam.gridss.working/chm13.bam.$metric_suffix $metrics_prefix.$metric_suffix
				fi
			done
			if [[ ! -f gen/$vcf ]] ; then
				cat > gen/slurm_virusbreakend_${n}_${depth}x.sh << EOF
#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=12:00:00
#SBATCH -p regular
#SBATCH --output=log.virusbreakend.%x.out
#SBATCH --error=log.virusbreakend.%x.err
. ~/conda_crest/etc/profile.d/conda.sh
conda activate virusbreakend
export PATH=~/projects/virusbreakend/gridss/scripts:\$PATH
export PATH=~/projects/virusbreakend/gridss/src/main/c/gridsstools:\$PATH
cd $PWD
cd gen
virusbreakend.sh --db $virusbreakenddb -r $ref -t 4 -o $vcf -j $GRIDSS_JAR --rmargs "-e rmblast" -w $PWD/gen/ $bam
EOF
			else
				echo "Found $vcf - skipping virusbreakend"
			fi
		fi
		mkdir -p gen/batvi_${n}_${depth}x
		if [[ ! -f gen/batvi_${n}_${depth}x/final_hits.txt ]] ; then
			if [[ ! -f gen/batvi_${n}_${depth}x/filelist.txt ]] ; then
				echo "$file.${depth}x.1.fq;$file.${depth}x.1.fq;500" > gen/batvi_${n}_${depth}x/filelist.txt
				echo > gen/batvi_${n}_${depth}x/batvirun.log.placeholder
			fi
			cat > gen/slurm_batvi_${n}_${depth}x.sh << EOF
#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=12:00:00
#SBATCH -p regular
#SBATCH --output=log.batvi.%x.out
#SBATCH --error=log.batvi.%x.err
. ~/conda_crest/etc/profile.d/conda.sh
conda activate virusbreakend
cd $PWD
cd batvi_${n}_${depth}x
echo > batvirun.log.placeholder
~/projects/virusbreakend/batvi/batvi1.03/call_integrations.sh $PWD/gen/batvi_${n}_${depth}x
EOF
		fi
		verse_dir=gen/verse/${n}/${depth}x
		if [[ ! -f $verse_dir/integration-sites.txt ]] ; then
			if [[ ! -d $verse_dir ]] ; then
				mkdir -p $verse_dir
				cp ~/projects/virusbreakend/verse/verse.config $verse_dir/verse.config
				echo "alignment_file = $bam" >> $verse_dir/verse.config
			fi
			cat > gen/slurm_verse_${n}_${depth}x.sh << EOF
#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=10g
#SBATCH --time=48:00:00
#SBATCH -p regular
#SBATCH --output=$(realpath $verse_dir)/log.verse.%x.out
#SBATCH --error=$(realpath $verse_dir)/log.verse.%x.err
. ~/conda_crest/etc/profile.d/conda.sh
conda activate verse
cd $(realpath $verse_dir)
perl ~/projects/virusbreakend/verse/VirusFinder2.0/VirusFinder.pl -c verse.config
EOF
		fi
		vifi_dir=gen/vifi/${n}/${depth}x
		mkdir -p $vifi_dir
		export VIFI_DIR=~/projects/virusbreakend/ViFi_hbv_hg19
		export AA_DATA_REPO=$VIFI_DIR/data_repo
		export REFERENCE_REPO=$VIFI_DIR/data
		if [[ ! -f $vifi_dir/output.clusters.txt.range ]] ; then
			cat > gen/slurm_vifi_${n}_${depth}x.sh << EOF
#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=10g
#SBATCH --time=48:00:00
#SBATCH -p regular
#SBATCH --output=$(realpath $vifi_dir)/log.vifi.%x.out
#SBATCH --error=$(realpath $vifi_dir)/log.vifi.%x.err
. ~/conda_crest/etc/profile.d/conda.sh
conda activate vifi
export VIFI_DIR=$VIFI_DIR
export AA_DATA_REPO=$AA_DATA_REPO
export REFERENCE_REPO=$REFERENCE_REPO
cd $(realpath $vifi_dir)
python $VIFI_DIR/scripts/run_vifi.py -f $file.${depth}x.1.fq -r $file.${depth}x.2.fq -v hbv -c 4
#-l $REFERENCE_REPO/hbv/hmms/hmms.txt
# -b $bam
EOF
		fi
		gridss_dir=gen/gridss/${n}/${depth}x
		mkdir -p $gridss_dir
		if [[ ! -f $gridss_dir/ ]] ; then
			cat > gen/slurm_gridss_${n}_${depth}x.sh << EOF
#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --time=48:00:00
#SBATCH -p regular
#SBATCH --output=$(realpath $vifi_dir)/log.gridss.%x.out
#SBATCH --error=$(realpath $vifi_dir)/log.gridss.%x.err
. ~/conda_crest/etc/profile.d/conda.sh
conda activate virusbreakend

# TODO: 
cd $(realpath $vifi_dir)
copy WGS metrics
gridss.bam -r $ref -a assembly.bam -o gridss.raw.vcf $bam
kraken annotate
filter to viral insertions
EOF
		fi
	done
done
chmod +x gen/*.sh

# merge outputs
#cat $(find . -name 'final_hits.txt') | grep -v LIB > gen/all_final_hits.txt
#grep MSA gen/all_final_hits.txt > ../publicdata/sim/batvi_all_final_hits_MSA.txt
#grep -v MSA gen/all_final_hits.txt > ../publicdata/sim/batvi_all_final_hits_noMSA.txt
#grep "^#" $(find gen/gen -name 'virusbreakend*.vcf' | head -1) > ../publicdata/sim/virusbreakend.vcf
#grep -v "##" $(find gen/gen -name 'virusbreakend*.vcf') | grep -v "#CHROM" >> ../publicdata/sim/virusbreakend.vcf
#grep -v Confidence $(find gen/verse -name integration-sites.txt) > ../publicdata/sim/verse.tsv
#grep chr $(find gen/vifi -name 'output.clusters.txt.range') | tr ':' ',' > ../publicdata/sim/vifi.tsv
