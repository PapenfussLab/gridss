#!/bin/bash
#
# runs breakdancer against bams
#
. common.sh
DINDEL=$BASE_DIR/tools/dindel/binaries/dindel-1.01-linux-64bit
DINDELPY_DIR=$BASE_DIR/tools/dindel/dindel-1.01-python
PATH=$DINDELPY_DIR:$PATH

# TODO: why no output?  # http://seqanswers.com/forums/showthread.php?t=17131
for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	CX_BAM=$BAM
	CX_CALLER=dindel
	cx_save
	ID=$(basename CX)
	XC_OUTPUT=$CX.vcf
	# single-threaded & takes ages: try on non-cluster nodes
	XC_NOCLUSTER=1
	XC_SCRIPT="rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	$DINDEL --analysis getCIGARindels --bamFile $BAM --outputFile $CX/$ID.dindel_output --ref $CX_REFERENCE
	makeWindows.py --inputVarFile $CX/$ID.dindel_output.variants.txt --windowFilePrefix $CX/$ID.realign_windows --numWindowsPerFile 1000
	for WINDOW_FILE in $CX/$ID.realign_windows.*.txt ; do
		$DINDEL --analysis indels --doDiploid --bamFile $BAM --ref $CX_REFERENCE --varFile \$WINDOW_FILE --libFile $CX/$ID.dindel_output.libraries.txt --outputFile \$WINDOW_FILE.stage2
	done
	ls -1 *.glf.txt | mergeOutputDiploid.py --inputFiles /dev/stdin --outputFile $CX.vcf --ref $CX_REFERENCE
	"
	xc_exec
done

