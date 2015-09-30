#!/bin/bash
# calculates the rate of compound assembly


# $1 bam suffix and output filename
calculate() {
	UNMAPPED=$( samtools merge -c -p -f /dev/stdout *.realign.$1.bam | samtools view -f 4 -c - )
	MAPPED=$(samtools merge -c -p -f /dev/stdout *.realign.$1.bam | samtools view -F 4 -c - )
	echo $1	$MAPPED	$UNMAPPED
}
echo RealignIteration	Mapped	Unmapped
calculate 0
calculate 1
calculate 2

