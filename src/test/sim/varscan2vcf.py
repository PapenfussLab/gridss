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
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""

i = 0
for line in sys.stdin:
	i += 1
	if i == 1: continue #skip header
	#Chrom   Position        Ref     Var     Cons:Cov:Reads1:Reads2:Freq:P-value     StrandFilter:R1+:R1-:R2+:R2-:pval       SamplesRef      SamplesHet      SamplesHom      SamplesNC       Cons:Cov:Reads1:Reads2:Freq:P-value     no_readcounts
	input = line.split('	')
	chr = input[0]
	pos = int(input[1])
	ref = input[2]
	alt = input[3]
	cons = input[4]
	strandfilter = input[5]
	# ignore the rest of the output for now
	svLen = len(alt) - 1
	if alt[0] == "+":
		svType = "INS"
		alt = ref + alt[1:]
		end = pos + svLen
	else:
		svType = "DEL"
		svLen *= -1
		alt = alt[1:]
		ref = ref + alt
		end = pos

	print "{0}	{1}	varscan{2}	{3}	{4}	{5}	PASS	END={6};SVLEN={7};SVTYPE={8}".format(
		chr,
		pos,
		i,
		ref,
		alt,
		".", # can we calcualte QUAL from P-value directly?
		end,
		svLen,
		svType
		)
	


