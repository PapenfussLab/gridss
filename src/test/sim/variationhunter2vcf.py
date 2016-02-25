#!/usr/bin/env python
#
# converts variationhunter output into VCF
#
import sys

print """##fileformat=VCFv4.1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=CILEN,Number=1,Type=Integer,Description="">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="">
##INFO=<ID=CILEN,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##ALT=<ID=DEL,Description="">
##ALT=<ID=INS,Description="">
##ALT=<ID=INV,Description="">
##INFO=<ID=chr,Number=1,Type=String,Description="The chromosome">
##INFO=<ID=Start_Outer,Number=1,Type=String,Description="the outer left coordinate of the mappings">
##INFO=<ID=Start_Inner,Number=1,Type=String,Description="the inner left coordinate of the mappings">
##INFO=<ID=End_Inner,Number=1,Type=String,Description="the inner right coordinate of the mappings">
##INFO=<ID=End_Outer,Number=1,Type=String,Description="the outer right coordinate of the mappings">
##INFO=<ID=sup,Number=1,Type=String,Description="total number of reads supporting the SV in all the individuals">
##INFO=<ID=VH_SVtype,Number=1,Type=String,Description="The type of SV (D : Deletion; V: inVersion; and I : Insertions)">
##INFO=<ID=AvgEditDits,Number=1,Type=String,Description="Average edit distance of all the mappings supporting the SV">
##INFO=<ID=MobileName,Number=1,Type=String,Description="">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""

def lineToDict(line):
	splitline = line.split(" ")
	return dict(zip(map(lambda s: s.split(":")[0], splitline), map(lambda s: s.split(":")[1], splitline)))

def process(line, i):
	ld = lineToDict(line)
	Start_Outer = int(ld["Start_Outer"])
	Start_Inner = int(ld["Start_Inner"])
	End_Inner = int(ld["End_Inner"])
	End_Outer = int(ld["End_Outer"])
	pos = (Start_Outer + Start_Inner) / 2
	endpos = (End_Outer + End_Inner) / 2
	if ld["SVtype"] == "D" : # Deletion
		len = pos - endpos # VCF SVLEN for DEL is negative
		out = "{0}	{1}	line{2}	N	<DEL>	.	PASS	SVTYPE=DEL;SVLEN={3};CILEN={4},{5};CIPOS={6},{7};CIEND={8},{9}".format(
			ld["Chr"], pos, i,
			len,
			len + int(ld["maxDelLen"]), len + int(ld["minDelLen"]), 
			Start_Outer - pos, Start_Inner - pos,
			End_Inner - endpos, End_Outer - endpos
			)
	elif ld["SVtype"] == "V" : # Inversion
		out = "TODO"
	elif ld["SVtype"] == "I" : #Insertion
		out = "TODO"
	elif ld["SVtype"] == "A" :
		# //For Mobile Insertions (chrN) the startPos is for reference genome and stopPos is for chrN
		# // Insertion from chrN in Forward orientation
		out = "TODO"
	elif ld["SVtype"] == "B" :
		# // Insertion from chrN in Reverse orientation
		out = "TODO"
	else:
		# ERROR - unknown variant type
		out = "unknown variant type " + ld["SVtype"]
	# Common suffixies
	print out + ";IMPRECISE;sup={0};VH_SVtype={1};AvgEditDits={2}".format(ld["sup"], ld["SVtype"], ld["AvgEditDits"])
	
i = 0
for line in sys.stdin:
	i = i + 1
	if line.startswith("Chr:") :
		process(line, i)



