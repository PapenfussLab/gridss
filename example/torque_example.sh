#!/bin/bash
base_dir=/path/to/project/
ref=$base_dir/ref/reference_genome.fa
working_dir=$base_dir/working
inputs="$base_dir/input1.bam $base_dir/input2.bam $base_dir/input3.bam"
inputarr=($inputs)
sample_name=example_torque
job_prefix=gridss_${sample_name}

mkdir -p $working_dir
cat > $working_dir/gridss.properties << EOF
# The default chunk size of 10Mb only gives ~300 units of work
# reducing to 1Mb gives ~3000 which gives a more even distribution
# of work across jobs
chunkSize=1000000
chunkSequenceChangePenalty=100000
EOF

#
# WARNING: assumes GRIDSS setupreference has already been run
#

# Assumes gridss.sh release artifacts are in $base_dir
gridss_cmd_common="$base_dir/gridss.sh \
	-r $ref \
	-o $base_dir/$sample_name.sv.vcf \
	-a $base_dir/$sample_name.asm.bam \
	-j $base_dir/gridss-2.6.3-gridss-jar-with-dependencies.jar \
	-w $working_dir \
	-c $working_dir/gridss.properties"

# Create a job array for each input file
threads=8
cat > $working_dir/${job_prefix}_preprocess.sh << EOF
#!/bin/bash
#PBS -N ${job_prefix}_preprocess
##PBS -o ${job_prefix}_preprocess
#PBS -t 1-${#inputarr[@]}
#PBS -l nodes=1:ppn=$threads,mem=16gb,walltime=144:00:0
#PBS -j oe

module add bwa samtools R java
inputarr=($inputs)
cd $working_dir
$gridss_cmd_common -t $threads \
	-s preprocess \
	\${inputarr[\$((PBS_ARRAYID -1 ))]}
EOF

# Create a job array containing assembly jobs
assembly_jobs=32
cat > $working_dir/${job_prefix}_assembly.sh << EOF
#!/bin/bash
#PBS -N ${job_prefix}_assembly
##PBS -o ${job_prefix}_assembly
#PBS -t 1-$assembly_jobs
#PBS -l nodes=1:ppn=$threads,mem=30gb,walltime=144:00:0
#PBS -j oe
module add bwa samtools R java
inputarr=($inputs)
cd $working_dir
$gridss_cmd_common -t $threads \
	-s assemble \
	--jobindex \$((PBS_ARRAYID -1 )) \
	--jobnodes $assembly_jobs \
	$inputs
EOF

# And finally, finish perform assembly gather
# variant calling, and annotation
threads=16
cat > $working_dir/${job_prefix}_call.sh << EOF
#!/bin/bash
#PBS -N ${job_prefix}_call
##PBS -o ${job_prefix}_call
#PBS -l nodes=1:ppn=$threads,mem=30gb,walltime=144:00:0
#PBS -j oe
module add bwa samtools R java
inputarr=($inputs)
cd $working_dir
$gridss_cmd_common -t $threads \
	$inputs
EOF

# queue the jobs up with a job dependency on the previous step
echo "pbs_id=\$(qsub $working_dir/${job_prefix}_preprocess.sh) ; echo \$pbs_id"
echo "pbs_id=\$(qsub -W depend=afterokarray:\$pbs_id $working_dir/${job_prefix}_assembly.sh) ; echo \$pbs_id)"
echo "qsub -W depend=afterokarray:\$pbs_id $working_dir/${job_prefix}_call.sh"

