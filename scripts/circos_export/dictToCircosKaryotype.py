#!/usr/bin/python
#
# converts gridss BEDPE output to circos links
#
import sys
import os

i = 0
for line in sys.stdin:
	i += 1
	if not line.startswith("@SQ"):
		continue
	input = line.split('	')
	contig = input[1].lstrip("SN:")
	length = int(input[2].lstrip("LN:"))
	
	# Use UCSC color scheme
	color = contig.lstrip("chr").split("_")[0]
	color = "chr" + color
	
	print "chr - {0} {1} {2} {3} {4}".format(contig, contig.lstrip("chr"), 0, length, color)
