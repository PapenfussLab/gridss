# Publishes the GRIDSS docker
unzip gridss_minimal_reproduction_data_for_error_1.zip
samtools faidx masked_ref.fa &
bwa index masked_ref.fa &
for f in *.gridss.working ; do
	bam=${f/.gridss.working/}
	cp $f/$bam.sv.bam $bam
	samtools index $bam &
done
wait