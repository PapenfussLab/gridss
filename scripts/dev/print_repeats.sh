#!/bin/bash
sambamba sort -q -n $1 -o $1.print_repeats.name_sorted.bam
samtools view $1.print_repeats.name_sorted.bam | uniq -d | tee $1.print_repeats.duplicates.sam