#!/bin/bash
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
	wgsim -r 0 -R 0 -X 0 -N 20000000 -1 $readlen -2 $readlen chm13.draft_v1.0.fasta gen/chm13.1.fq gen/chm13.2.fq
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
			wgsim -r 0 -R 0 -X 0 -N $(( $depth * $reflen / $pairbases )) -1 150 -2 150 $file $file.${depth}x.1.fq $file.${depth}x.2.fq
		fi
		vcf=gen/virusbreakend_${n}_${depth}x.vcf
		if [[ ! -f $vcf.virusbreakend.working/$(basename $vcf).kraken2.report.viral.extracted.txt ]] ; then
			bam=$file.${depth}x.bam
			metrics_prefix=$vcf.virusbreakend.working/adjusted/$(basename $bam).viral.bam.gridss.working/$(basename $bam).viral.bam
			mkdir -p $(dirname $metrics_prefix)
			for metric_suffix in cigar_metrics idsv_metrics insert_size_metrics mapq_metrics tag_metrics ; do
				if [[ ! -f $metrics_prefix.$metric_suffix ]] ; then
					cp gen/chm13.bam.gridss.working/chm13.bam.$metric_suffix $metrics_prefix.$metric_suffix
				fi
			done
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
if [[ ! -f $bam ]] ; then
	echo "Creating $bam"
	bwa mem -t 4 $ref $file.${depth}x.1.fq $file.${depth}x.2.fq | samtools sort -@ 4 -o $bam -O BAM -
	samtools index $bam
fi
virusbreakend.sh --db $virusbreakenddb -r $ref -t 4 -o $vcf -j $GRIDSS_JAR --rmargs "-e rmblast" -w $PWD/gen/ $bam
EOF
		fi
		mkdir -p gen/batvi_${n}_${depth}x
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
		verse_dir=gen/verse/${n}/${depth}x
		if [[ ! -d $verse_dir ]] ; then
			mkdir -p $verse_dir
			cp ~/projects/virusbreakend/verse/verse.config $verse_dir/verse.config
			echo "alignment_file = $bam" >> $verse_dir/verse.config
		fi
		cat gen/slurm_verse_${n}_${depth}x.sh << EOF
#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=48:00:00
#SBATCH -p regular
#SBATCH --output=$verse_dir/log.verse.%x.out
#SBATCH --error=$verse_dir/log.verse.%x.err
. ~/conda_crest/etc/profile.d/conda.sh
conda activate virusbreakend
cd $verse_dir
perl ~/projects/virusbreakend/verse/VirusFinder2.0/VirusFinder.pl -c $verse_dir/verse.config -o $verse_dir
EOF
	done
done
chmod +x gen/*.sh

# merge outputs
#cat $(find . -name 'final_hits.txt') | grep -v LIB > gen/all_final_hits.txt
#grep MSA gen/all_final_hits.txt > ../publicdata/sim/batvi_all_final_hits_MSA.txt
#grep -v MSA gen/all_final_hits.txt > ../publicdata/sim/batvi_all_final_hits_noMSA.txt
#grep "^#" $(find gen/gen -name 'virusbreakend*.vcf' | head -1) > ../publicdata/sim/virusbreakend.vcf
#grep -v "##" $(find gen/gen -name 'virusbreakend*.vcf') | grep -v "#CHROM" >> ../publicdata/sim/virusbreakend.vcf

