#!/usr/bin/env python
#
# converts varscan v2.3.6 output into valid VCF
#
import sys

# TODO: custom INFO fields for varscan2-specific values
print """##fileformat=VCFv4.1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=CIPOS,Number=2,Type=Integer>
##INFO=<ID=CIRPOS,Number=2,Type=Integer>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""

i = 0
for line in sys.stdin:
	i += 1
	#Chrom   Position        Ref     Var     Cons:Cov:Reads1:Reads2:Freq:P-value     StrandFilter:R1+:R1-:R2+:R2-:pval       SamplesRef      SamplesHet      SamplesHom      SamplesNC       Cons:Cov:Reads1:Reads2:Freq:P-value     no_readcounts
	input = line.split('	')
	chr1 = input[0]
	start1 = int(input[1])
	end1 = int(input[2])
	chr1 = input[3]
	start2 = int(input[4])
	end2 = int(input[5])
	id = input[6]
	#size = int(input[7])
	strand1 = input[8]
	strand2 = input[9]
	type = "TYPE:DELETION" #input[10]
	if type.startswith("TYPE:DELETION") :
		print "{0}	{1}	{2}	{3}	{4}	{5}	PASS	END={6};SVLEN={7};SVTYPE=DEL;CIPOS={8},{9};CIRPOS={10},{11}".format(
			chr1,
			start1,
			id,
			"N",
			"<DEL>",
			".",
			start2,
			end2 - end1, #size,
			0, end1 - start1,
			0, end2 - start2
			)
	


