#!/usr/bin/env Rscript
#
# Sanity checks a GRIDSS VCF file to ensure that the output makes sense
#
usage <- "Usage: gridsssanitycheck <input VCF>"
args <- commandArgs(TRUE)
if (length(args) != 2) {
	write(usage, stderr())
	q(save="no", status=1)
}
if (!file.exists(args[2])) {
	write(paste(args[2], "not found"), stderr())
	q(save="no", status=1)
}
# quieten all the dependency package conflict spam
library(BiocGenerics, quietly=TRUE, warn.conflicts=FALSE)
library(S4Vectors, quietly=TRUE, warn.conflicts=FALSE)
library(IRanges, quietly=TRUE, warn.conflicts=FALSE)
library(matrixStats, quietly=TRUE, warn.conflicts=FALSE)
library(DelayedArray, quietly=TRUE, warn.conflicts=FALSE)
library(XVector, quietly=TRUE, warn.conflicts=FALSE)
library(Biostrings, quietly=TRUE, warn.conflicts=FALSE)
suppressPackageStartupMessages(library(Biobase))
# Ok, now load direct package dependencies. It's a pity warn.conflicts does not get passed through
library(VariantAnnotation, quietly=TRUE, warn.conflicts=FALSE)
#library(rtracklayer, quietly=TRUE, warn.conflicts=FALSE)
#library(tidyverse, quietly=TRUE, warn.conflicts=FALSE)
#library(stringr, quietly=TRUE, warn.conflicts=FALSE)
#library(testthat, quietly=TRUE, warn.conflicts=FALSE)
#library(stringdist, quietly=TRUE, warn.conflicts=FALSE)

errorCount <- 0
margin <- 0.04 # margin of error allowed due to VCF rounding to 2dp
vcf <- readVcf(args[2], "unknown")

vcf = vcf[elementNROWS(info(vcf)$MATEID) == 1] # TODO expand script to include single breakends

######
# Annotation breakdown validations

# QUAL INFO breakdown
mismatchedField = abs(rowRanges(vcf)$QUAL - (info(vcf)$ASQ + info(vcf)$RASQ + info(vcf)$CASQ + info(vcf)$SRQ + info(vcf)$IQ + info(vcf)$RPQ)) > margin
if (any(mismatchedField)) {
	errorCount <- errorCount + 1
	write(paste0("INFO quality score breakdown does not match total for: ", paste(row.names(vcf)[mismatchedField], collapse=", ")), stdout())
}
# QUAL FORMAT breakdown
mismatchedField = rowSums(abs(geno(vcf)$ASQ + geno(vcf)$RASQ + geno(vcf)$CASQ + geno(vcf)$SRQ + geno(vcf)$IQ + geno(vcf)$RPQ - geno(vcf)$QUAL)) > margin
if (any(mismatchedField)) {
	errorCount <- errorCount + 1
	write(paste0("FORMAT quality score breakdown does not match FORMAT QUAL total for: ", paste(row.names(vcf)[mismatchedField], collapse=", ")), stdout())
}
# BQ INFO breakdown
mismatchedField = abs(info(vcf)$BQ - (info(vcf)$BAQ + info(vcf)$BSCQ + info(vcf)$BUMQ)) > margin
if (any(mismatchedField)) {
	errorCount <- errorCount + 1
	write(paste0("INFO BQ breakend quality score breakdown does not match total for: ", paste(row.names(vcf)[mismatchedField], collapse=", ")), stdout())
}
# BQ FORMAT breakdown
mismatchedField = rowSums(abs(geno(vcf)$BAQ + geno(vcf)$BSCQ + geno(vcf)$BUMQ - geno(vcf)$BQ)) > margin
if (any(mismatchedField)) {
	errorCount <- errorCount + 1
	write(paste0("FORMAT quality score breakdown does not match FORMAT QUAL total for: ", paste(row.names(vcf)[mismatchedField], collapse=", ")), stdout())
}


# Breakends are partnered
missingMATEID <- row.names(vcf)[!(as.character(info(vcf)$MATEID) %in% row.names(vcf))]
if (length(missingMATEID) > 0 ) {
	errorCount <- errorCount + 1
	write(paste("Missing partner specified in MATEID for:", paste(missingMATEID, collapse=", ")), stdout())
	vcf <- vcf[!(row.names(vcf) %in% missingMATEID),]
}

pvcf <- vcf[as.character(info(vcf)$MATEID),]

mismatchedField = rowRanges(vcf)$QUAL != rowRanges(pvcf)$QUAL
if (any(mismatchedField)) {
	errorCount <- errorCount + 1
	write(paste0("QUAL does not match partner for: ", paste(row.names(vcf)[mismatchedField], collapse=", ")), stdout())
}
mismatchedField = rowRanges(vcf)$FILTER != rowRanges(pvcf)$FILTER
if (any(mismatchedField)) {
	errorCount <- errorCount + 1
	write(paste0("FILTER does not match partner for: ", paste(row.names(vcf)[mismatchedField], collapse=", ")), stdout())
}

breakpointFieldName <- c("CAS", "AS", "RAS", "SR", "RP", "IC", "ASRP", "ASSR", "ASQ", "RASQ", "SRQ", "RPQ", "IQ", "CQ")
breakpointFieldPartner <- c("CAS", "RAS", "AS", "SR", "RP", "IC", "ASRP", "ASSR", "RASQ", "ASQ", "SRQ", "RPQ", "IQ", "CQ")
# INFO matches partner
mismatchedFieldList = data.frame(field=character(0), count=numeric(0))
for (i in seq_along(breakpointFieldName)) {
	mismatchedField = abs(info(vcf)[[breakpointFieldName[i]]] - info(pvcf)[[breakpointFieldPartner[i]]]) > margin
	if (any(mismatchedField)) {
		errorCount <- errorCount + 1
		write(paste0(breakpointFieldName[i], ": INFO does not match partner for: ", paste(row.names(vcf)[mismatchedField], collapse=", ")), stdout())
		mismatchedFieldList = rbind(mismatchedFieldList, data.frame(field=breakpointFieldName[i], sum(mismatchedField)))
	}
}
# FORMAT matches partner
for (i in seq_along(breakpointFieldName)) {
	if (breakpointFieldName[i] %in% c("AS", "RAS", "CAS", "CQ")) {
		# Ignore AS and RAS as they're not output
		next()
	}
	mismatchedField = abs(rowSums(geno(vcf)[[breakpointFieldName[i]]] - geno(pvcf)[[breakpointFieldPartner[i]]])) > margin
	if (any(mismatchedField)) {
		errorCount <- errorCount + 1
		write(paste0(breakpointFieldName[i], ": FORMAT does not match partner for: ", paste(row.names(vcf)[mismatchedField], collapse=", ")), stdout())
	}
}
# FORMAT matches INFO
for (fieldName in c("SR", "RP", "ASRP", "ASSR", "RASQ", "ASQ", "SRQ", "RPQ", "REF", "REFPAIR", "BQ", "BUM", "BSC", "BAQ", "BUMQ", "BSCQ")) {
	mismatchedField = abs(info(vcf)[[fieldName]] - rowSums(geno(vcf)[[fieldName]])) > margin
	if (any(mismatchedField)) {
		errorCount <- errorCount + 1
		write(paste0(fieldName, ": INFO does not match FORMAT for: ", paste(row.names(vcf)[mismatchedField], collapse=", ")), stdout())
	}
}
# CIPOS matches CIRPOS
mismatchedField = all(info(vcf)$CIPOS != info(pvcf)$CIRPOS)
mismatchedField = !is.na(mismatchedField) & mismatchedField
if (any(mismatchedField)) {
	errorCount <- errorCount + 1
	write(paste0("CIPOS does not match partner CIRPOS for: ", paste(row.names(vcf)[mismatchedField], collapse=", ")), stdout())
}
mismatchedField = any(info(vcf)$HOMLEN - info(pvcf)$HOMLEN > 1)
if (any(mismatchedField)) {
	errorCount <- errorCount + 1
	write(paste0("IHOMPOS length does not match partner for: ", paste(row.names(vcf)[mismatchedField], collapse=", ")), stdout())
}

q(save="no", status=errorCount)
