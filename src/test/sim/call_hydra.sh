 #!/bin/bash
#
# runs breakdancer against bams
#
. common.sh
export PATH=$BASE_DIR/tools/hydra-multi:$PATH
export PATH=$BASE_DIR/tools/hydra-multi/bin:$PATH
export PATH=$BASE_DIR/tools/hydra-multi/scripts:$PATH

# https://github.com/chapmanb/bcbb/blob/master/nextgen/scripts/hydra_to_vcf.py
HYDRA2VCF=$BASE_DIR/tools/hydra_to_vcf.py

for BAM in $DATA_DIR/*.sq.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=hydra/master20160129
	cx_save
	MAD=`echo "scale=4; $CX_READ_FRAGMENT_STDDEV / 1.4826" | bc`
	MLD=`echo "scale=4; 10 * $MAD" | bc`
	MNO=`echo "scale=4; 20 * $MAD" | bc`
	# follow the hydra workflow from https://code.google.com/p/hydra-sv/wiki/TypicalWorkflow
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="rm -rf $CX; mkdir $CX $CX 2>/dev/null; cd $CX
ulimit -n 16384
ulimit -n 50000
ln -s input.bam $CX_BAM
echo 'sample	input.bam' > config.stub.txt
hydra-multi.sh run -t $XC_CORES config.stub.txt
$BASE_DIR/hydra2vcf.py < $CX/all.hydra.sv.final > $CX.vcf
	"
	xc_exec
done

