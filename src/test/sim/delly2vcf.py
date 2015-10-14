#!/usr/local/bioinf/bin/python
#
# converts delly output into pseudo-VCF
# Delly format: http://www.embl.de/~rausch/delly.html

# TODO: use breakpoint concensus
# TODO: process invy output
# TODO: process jmpy output
import sys
import os

print """##fileformat=VCFv4.1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=NPE,Number=1,Type=Integer,Description="Number of supporting paired-end reads">
##INFO=<ID=NSR,Number=1,Type=Integer,Description="Number of supporting split reads">
##INFO=<ID=AMQ,Number=1,Type=Float,Description="Average mapping quality of pair-end reads">
##INFO=<ID=PCAQVR,Number=1,Type=Float,Description="Percent consensus alignment quality vs reference">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">
##ALT=<ID=CNV,Description="Copy number variable region">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""

#chr, start, end, size, #supporting_pairs, avg._mapping_quality, deletion_id
srSupport = {}
if os.path.exists(sys.argv[1] + ".del.br.txt"):
	for line in open( sys.argv[1] + ".del.br.txt", "r" ):
		if line.find("Deletion_") == -1: continue
		input = line.split('	')
		id = input[6].strip().strip("<>")
		srSupport[id] = input
for line in open( sys.argv[1] + ".del.txt", "r" ):
	if line.find("Deletion_") == -1: continue
	input = line.split('	')
	chr = input[0]
	start = int(input[1])
	end = int(input[2])
	size = int(input[3])
	peCount = int(input[4])
	avgMappingQual = float(input[5])
	id = input[6].strip().strip("<>")
	srCount = 0
	percentConsensusAlignmentQualityVsReference = 0
	if srSupport.has_key(id):
		input = srSupport[id]
		# use precise location
		start = int(input[1])
		end = int(input[2])
		size = int(input[3])
		# add new fields
		srCount = input[4]
		percentConsensusAlignmentQualityVsReference = input[5]
	print "{0}	{1}	{2}	N	<DEL>	{3}	PASS	END={4};SVLEN={5};SVTYPE={6};NPE={7};AMQ={8};NSR={9};PCAQVR={10}{11}".format(
			chr,
			start,
			id,
			0,
			end,
			-size,
			"DEL",
			peCount,
			avgMappingQual,
			srCount,
			percentConsensusAlignmentQualityVsReference,
			";IMPRECISE" if srCount == 0 else ""
			)
			
#chr, start, end, size, #supporting_pairs, avg._mapping_quality, tandem_duplication_id
srSupport = {}
if os.path.exists(sys.argv[1] + ".dup.br.txt"):
	for line in open( sys.argv[1] + ".dup.br.txt", "r" ):
		if line.find(">Duplication_") == -1: continue
		input = line.split('	')
		id = input[6].strip().strip("<>")
		srCount = input[4]
		percentConsensusAlignmentQualityVsReference = input[5]
		srSupport[id] = (srCount, percentConsensusAlignmentQualityVsReference)
for line in open( sys.argv[1] + ".dup.txt", "r" ):
	if line.find(">Duplication_") == -1: continue
	input = line.split('	')
	chr = input[0]
	start = int(input[1])
	end = int(input[2])
	size = int(input[3])
	peCount = int(input[4])
	avgMappingQual = float(input[5])
	id = input[6].strip().strip("<>")
	srCount = 0
	percentConsensusAlignmentQualityVsReference = 0
	if srSupport.has_key(id):
		srCount = srSupport[id][0]
		percentConsensusAlignmentQualityVsReference = srSupport[id][1]
	print "{0}	{1}	{2}	N	<DUP:TANDEM>	{3}	PASS	END={4};SVLEN={5};SVTYPE={6};NPE={7};AMQ={8};NSR={9};PCAQVR={10}{11}".format(
			chr,
			start,
			id,
			0,
			end,
			size,
			"DUP",
			peCount,
			avgMappingQual,
			srCount,
			percentConsensusAlignmentQualityVsReference,
			";IMPRECISE" if srCount == 0 else ""
			)

