#!/bin/bash
#conda create -n vifi bwa=0.7.17 python=2.7 pysam=0.15.2 samtools=1.9 hmmer
export PATH=~/dev/gridss/:$PATH
export GRIDSS_JAR=~/projects/virusbreakend/gridss-2.10.2-gridss-jar-with-dependencies.jar
hpv=MT783410.1
hbv=LC500247.1
chr1_len=248387497
n=0
pos=0
flank=50000
breakoffset=10
virus_length=2000
readlen=150
virusbreakenddb=~/projects/virusbreakend/virusbreakenddb
ref=~/projects/reference_genomes/human/hg19.fa
art_args="--noALN --paired --seqSys HSXn -ir 0 -ir2 0 -dr 0 -dr2 0 -k 0 -l 150 -m 500 -s 100 -rs 1 "
if [[ ! -f gen/chm13.bam.gridss.working/chm13.bam.insert_size_metrics ]] ; then
	file=$PWD/gen/chr13.fa
	depth=30
	if [[ ! -f $file ]] ; then
		echo ">chm13_chr1_LC500247" >> $file
		samtools faidx chm13.draft_v1.0.fasta chr1:1-33000000 | grep -v ">" | tr -d '\n' >> $file
		samtools faidx $hbv.fa $hbv:3215 | grep -v ">" | tr -d '\n' >> $file
		samtools faidx chm13.draft_v1.0.fasta chr1:33000001-248387497 | grep -v ">" | tr -d '\n' >> $file
		# GRIDSS uses WGS metrics which are going to be wrong for the small simulation - create a baseline
		art_illumina $art_args --fcov $depth --in $file  -o $file.${depth}x.
		cat ecoli_1.fq >> $file.${depth}x.1.fq
		cat ecoli_2.fq >> $file.${depth}x.2.fq
	fi
	if [[ ! -f gen/chm13.bam ]] ; then
		bwa mem -t $(nproc) $ref $file.${depth}x.1.fq $file.${depth}x.2.fq | samtools sort -@ $(nproc) -o gen/chm13.bam -O BAM -
	fi
	gridss.sh -r $ref -s preprocess gen/chm13.bam -w gen/
fi
exit 0
while [[ $n -lt 248 ]] ; do
	n=$(($n + 1))
	pos=$(($n * 1000000))
	virus_pos=$(($n * 4))
	mkdir -p $PWD/gen/${n}/
	file=$PWD/gen/${n}/chr1_${pos}.fa
	if [[ ! -f $file ]] ; then
		echo ">$file" >> $file
		samtools faidx chm13.draft_v1.0.fasta chr1:$(($pos - $flank))-$(($pos)) | grep -v ">" | tr -d '\n' >> $file
		samtools faidx $hbv.fa $hbv:$virus_pos-$(($virus_pos + $virus_length)) | grep -v ">" | tr -d '\n' >> $file
		samtools faidx chm13.draft_v1.0.fasta chr1:$(($pos + $breakoffset))-$(($pos + $flank + $breakoffset)) | grep -v ">" | tr -d '\n' >> $file
		echo "Generated $n"
	fi
	reflen=$(( 2 * $flank + $virus_length ))
	pairbases=$(( 2 * $readlen ))
	for depth in 5 10 15 30 60 ; do
		if [[ ! -f $file.${depth}x.1.fq ]] ; then
			conda activate virusbreakend
			art_illumina $art_args --fcov $depth --in $file -o $file.${depth}x.
			#wgsim -S 1 -r 0 -R 0 -X 0 -N $(( $depth * $reflen / $pairbases )) -1 150 -2 150 $file $file.${depth}x.1.fq $file.${depth}x.2.fq
			cat ecoli_1.fq >> $file.${depth}x.1.fq
			cat ecoli_2.fq >> $file.${depth}x.2.fq
		fi
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
echo "\$(date) Start"
. ~/conda_crest/etc/profile.d/conda.sh
conda activate virusbreakend
cd $PWD
cd gen
if [[ ! -f $bam ]] ; then
	echo "Creating $bam"
	bwa mem -t 4 $ref $file.${depth}x.1.fq $file.${depth}x.2.fq | samtools sort -@ 4 -o $bam -O BAM -
	samtools index $bam
fi
echo "\$(date) End"
EOF
		fi
		virusbreakend_dir=$PWD/gen/virusbreakend/${n}/${depth}x
		virusbreakend_vcf=$virusbreakend_dir/virusbreakend_${n}_${depth}x.vcf
		mkdir -p $virusbreakend_dir
		if [[ ! -f $virusbreakend_vcf.virusbreakend.working/$(basename $virusbreakend_vcf).kraken2.report.viral.extracted.txt ]] ; then
			metrics_prefix=$virusbreakend_vcf.virusbreakend.working/adjusted/$(basename $bam).viral.bam.gridss.working/$(basename $bam).viral.bam
			mkdir -p $(dirname $metrics_prefix)
			for metric_suffix in cigar_metrics idsv_metrics insert_size_metrics mapq_metrics tag_metrics ; do
				if [[ ! -f $metrics_prefix.$metric_suffix ]] ; then
					cp gen/chm13.bam.gridss.working/chm13.bam.$metric_suffix $metrics_prefix.$metric_suffix
				fi
			done
			if [[ ! -f $virusbreakend_vcf ]] ; then
				cat > gen/slurm_virusbreakend_${n}_${depth}x.sh << EOF
#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=12:00:00
#SBATCH -p regular
#SBATCH --output=log.virusbreakend.%x.out
#SBATCH --error=log.virusbreakend.%x.err
echo "\$(date) Start"
. ~/conda_crest/etc/profile.d/conda.sh
conda activate virusbreakend
export PATH=~/projects/virusbreakend/gridss/scripts:\$PATH
export PATH=~/projects/virusbreakend/gridss/src/main/c/gridsstools:\$PATH
cd $virusbreakend_dir
virusbreakend.sh --db $virusbreakenddb -r $ref -t 4 -o $virusbreakend_vcf -j $GRIDSS_JAR --rmargs "-e rmblast" -w . $bam
echo "\$(date) End"
EOF
			else
				echo "Found $virusbreakend_vcf - skipping virusbreakend"
			fi
		fi
		batvi_dir=$PWD/gen/batvi/${n}/${depth}x
		mkdir -p $batvi_dir
		if [[ ! -f $batvi_dir/final_hits.txt ]] ; then
			if [[ ! -f $batvi_dir/filelist.txt ]] ; then
				echo "$file.${depth}x.1.fq;$file.${depth}x.1.fq;500" > $batvi_dir/filelist.txt
				echo > $batvi_dir/batvirun.log.placeholder
			fi
			cat > gen/slurm_batvi_${n}_${depth}x.sh << EOF
#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=12:00:00
#SBATCH -p regular
#SBATCH --output=log.batvi.%x.out
#SBATCH --error=log.batvi.%x.err
echo "\$(date) Start"
. ~/conda_crest/etc/profile.d/conda.sh
conda activate virusbreakend
cd $batvi_dir
echo > batvirun.log.placeholder
~/projects/virusbreakend/batvi/batvi1.03/call_integrations.sh $batvi_dir
echo "\$(date) End"
EOF
		fi
		verse_dir=$PWD/gen/verse/${n}/${depth}x
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
echo "\$(date) Start"
. ~/conda_crest/etc/profile.d/conda.sh
conda activate verse
cd $(realpath $verse_dir)
perl ~/projects/virusbreakend/verse/VirusFinder2.0/VirusFinder.pl -c verse.config
echo "\$(date) End"
EOF
		fi
		vifi_dir=$PWD/gen/vifi/${n}/${depth}x
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
echo "\$(date) Start"
. ~/conda_crest/etc/profile.d/conda.sh
conda activate vifi
export VIFI_DIR=$VIFI_DIR
export AA_DATA_REPO=$AA_DATA_REPO
export REFERENCE_REPO=$REFERENCE_REPO
cd $(realpath $vifi_dir)
python $VIFI_DIR/scripts/run_vifi.py -f $file.${depth}x.1.fq -r $file.${depth}x.2.fq -v hbv -c 4
#-l $REFERENCE_REPO/hbv/hmms/hmms.txt
# -b $bam
echo "\$(date) End"
EOF
		fi
		gridss_dir=$PWD/gen/gridss/${n}/${depth}x
		gridss_vcf=$gridss_dir/gridss_${n}_${depth}x.vcf
		mkdir -p $gridss_dir
		if [[ ! -f $gridss_dir/ ]] ; then
			metrics_prefix=$gridss_dir/$(basename $bam).gridss.working/$(basename $bam)
			mkdir -p $(dirname $metrics_prefix)
			for metric_suffix in cigar_metrics idsv_metrics insert_size_metrics mapq_metrics tag_metrics ; do
				if [[ ! -f $metrics_prefix.$metric_suffix ]] ; then
					cp gen/chm13.bam.gridss.working/chm13.bam.$metric_suffix $metrics_prefix.$metric_suffix
				fi
			done
			if [[ ! -f $gridss_dir/$gridss_vcf ]]
			cat > gen/slurm_gridss_${n}_${depth}x.sh << EOF
#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=48:00:00
#SBATCH -p regular
#SBATCH --output=$(realpath $gridss_dir)/log.gridss.%x.out
#SBATCH --error=$(realpath $gridss_dir)/log.gridss.%x.err
echo "\$(date) Start"
. ~/conda_crest/etc/profile.d/conda.sh
conda activate virusbreakend
cd $(realpath $gridss_dir)
gridss.sh --jvmheap 14g -t 4 -j $GRIDSS_JAR -r $ref -a assembly.bam -o gridss.raw.vcf $bam
gridss_annotate_vcf_kraken2.sh \
	-o gridss.ann.vcf \
	-t 4 \
	--kraken2db $virusbreakenddb \
	gridss.raw.vcf
grep INSTAXID gridss.ann.vcf | grep -v "INSTAXID=9606;" > $gridss_vcf
echo "\$(date) End"
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
#grep "^#" gen/gridss/1/60x/gridss_1_60x.vcf > ../publicdata/sim/virusbreakend.vcf
#grep -v "##" $(find gen/gridss -name 'gridss_*x.vcf') | grep -v "#CHROM" | grep INSTAXID >> ../publicdata/sim/gridss2.vcf

