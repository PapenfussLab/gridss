#!/usr/bin/python
#
# converts gridss BEDPE output to circos links
#
import sys
import os
import math

i = 0
for line in sys.stdin:
	i += 1
	if line[0] == '#':
		# skip header line if exists
		continue
	input = line.split('	')
	chrom1 = input[0]
	start1 = int(input[1])
	end1 = int(input[2])
	chrom2 = input[3]
	start2 = int(input[4])
	end2 = int(input[5])
	name = input[6]
	score = float(input[7])
	strand1 = input[8]
	strand2 = input[9]
	hom = input[10]
	AS = int(input[11])
	RAS = int(input[12])
	SC = int(input[13])
	RSC = int(input[14])
	RP = int(input[15])
	ASQ = float(input[16])
	RASQ = float(input[17])
	SCQ = float(input[18])
	RSCQ = float(input[19])
	RPQ = float(input[20])
	
	color = "lgrey"
	if RAS + AS > 0:
		color="dgrey"
	if RAS > 0 and AS > 0:
		color="red"
	if score < 1000:
		color = color + "_a" + str(int(math.ceil(min(score, 1000) * 15 / 1000)) + 1)
	
	thickness=math.log10(score)
	if score >= 250:
		print "{0} {1} {2} {3} {4} {5} color={6},thickness={7}".format(
			chrom1, start1, end1, chrom2, start2, end2,
			color, thickness
			)
