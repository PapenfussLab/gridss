library(VariantAnnotation)
library(plyr)
source("../../test/r/libvcf.R")

gridss.annotateBreakpointHits <- function(bed, bedMate, gridssVcf, ...) {
  return (gridss.annotateBreakpoints(bed, bedMate, gridssVcf, ...)$bed)
}
gridss.annotateBreakpoints <- function(bedWithMate, gridssVcf, ...) {
  return (gridss.annotateBreakpoints(bedWithMate, bedWithMate[bedWithMate$mate,], gridssVcf))
}
gridss.annotateBreakpoints <- function(bed, bedMate, gridssVcf, ...) {
  # filter breakends and one-sided breakpoint calls
  gridssVcf <- gridssVcf[as.character(info(gridssVcf)$MATEID) %in% row.names(gridssVcf),]
  gridssdf <- gridss.vcftodf(gridssVcf)
  callPos <- rowRanges(gridssVcf)
  callPos$mate <- as.character(info(gridssVcf)$MATEID)
  strand(callPos) <- ifelse(str_detect(as.character(callPos$ALT), "[[:alpha:]]+(\\[|]).*(\\[|])"), "+", "-")
  callPosMate <- callPos[callPos$mate,]
  
  hits <- rbind(
    as.data.frame(findOverlaps(bed, callPos, ...)),
    as.data.frame(findOverlaps(bed, callPosMate, ...)),
    as.data.frame(findOverlaps(bedMate, callPos, ...)),
    as.data.frame(findOverlaps(bedMate, callPosMate, ...)))
  hits <- hits[duplicated(hits),]
  hits$qual <- rowRanges(gridssVcf)$QUAL[hits[[2]]]
  hits$hc <- gridssdf$QUAL[hits[[2]]] >= 1000 & gridssdf$AS[hits[[2]]] > 0 & gridssdf$RAS[hits[[2]]] > 0
  hits$gridssid <- names(rowRanges(gridssVcf))[hits[[2]]]
  hits$bedid <- names(bed)[hits[[1]]]
  hits <- hits[order(hits$qual),]
  
  annotatedBed <- bed
  annotatedBed$gridssid <- NA
  annotatedBed$gridssid[hits[[1]]] <- hits$gridssid
  annotatedBed$qual <- 0
  annotatedBed$qual[hits[[1]]] <- hits$qual
  annotatedGridss <- gridssdf
  annotatedGridss$bedid <- NA
  annotatedGridss$bedid[hits[[2]]] <- hits$bedid
  return(list(bed=annotatedBed, gridss=annotatedGridss))
}

gridss.overlaps <- function(vcf, bed, ...) {
  hits <- findOverlaps(rowRanges(vcf), bed, ...)
  result <- rep(FALSE, nrow(vcf))
  result[queryHits(hits)] <- TRUE
  return(result)
}
#CompressedIntegerList to array
ciltoarray <- function(cil, column) {
  if (length(cil) == 0) return(NULL)
  result <- matrix(unlist(cil), nrow=length(cil), byrow=TRUE)[,column]
  result[is.na(result)] <- 0
  return(result)
}
#CompressedList to matrix
cltomatrix <- function(cl) {
  el <- elementLengths(cl)
  rowCount <- length(cl)
  colCount <- max(el)
  # TODO: can we do with without having to use lapply?
  # colOffset <- rep(seq_len(rowCount), elementLengths(cl))
  # rowOffset <- #cumsum resetting
  mat <- matrix(unlist(lapply(cl, function(x) c(x, rep(0, colCount-length(x))))), rowCount, colCount, byrow=TRUE)
  return(mat)
}
# Flattens a list column into summary and index columns
gridss.vcftodf.flattenNumeric <- function(df, name, column, max=FALSE, allColumns=FALSE) {
  mat <- cltomatrix(column)
  
  if (max) {
    df[[name]] <- rowMax(mat)
  } else {
    df[[name]] <- rowSums(mat)
  }
  if (allColumns) {
    for (offset in seq_len(ncol(mat))) {
      df[paste0(name, offset)] <- mat[,offset]
    }
  }
  return(df)
}
gridss.removeUnpartnerededBreakend <- function(vcf) {
  toKeep <- as.character(info(vcf)$MATEID) %in% row.names(vcf)
  if (sum(toKeep == FALSE) > 0) {
    warning(paste("Found ", sum(toKeep == FALSE), " unpartnered breakends."))
  }
  return(vcf[toKeep,])
}
gridss.vcftodf <- function(vcf, allColumns=FALSE, sanityCheck=TRUE) {
  i <- info(vcf)
  df <- data.frame(variantid=names(rowRanges(vcf)))
  df$POS <- paste0(seqnames(rowRanges(vcf)), ":", start(rowRanges(vcf)))
  df$FILTER <- rowRanges(vcf)$FILTER
  df$QUAL <- rowRanges(vcf)$QUAL
  df$EVENT <-i$EVENT
  df$mate <- as.character(i$MATEID)
  df$SOMATIC <- i$SOMATIC
  matchLength <- as.integer(gsub(".*_", "", as.character(i$TRUTH_MATCHES)))
  mismatchLength <- as.integer(gsub(".*_", "", as.character(i$TRUTH_MISREALIGN)))
  df$SVLEN <- ifelse(is.na(matchLength), mismatchLength, matchLength)
  df$SVTYPE <- i$SVTYPE
  df$HOMSEQ <- ifelse(is.na(as.character(i$HOMSEQ)), "", as.character(i$HOMSEQ))
  df$HOMLEN <- as.numeric(i$HOMLEN)
  df$call <- ifelse(!is.na(matchLength), "good", ifelse(!is.na(mismatchLength), "misaligned", "bad"))
  df <- gridss.vcftodf.flattenNumeric(df, "REF", i$REF, allColumns=allColumns)
  df <- gridss.vcftodf.flattenNumeric(df, "REFPAIR", i$REFPAIR, allColumns=allColumns)
  df$SPV <- i$SPV
  df$CQ <- i$CQ
  df$BQ <- i$BQ
  
  df$AS <- i$AS
  df <- gridss.vcftodf.flattenNumeric(df, "RP", i$RP, allColumns=allColumns)
  df <- gridss.vcftodf.flattenNumeric(df, "SC", i$SC, allColumns=allColumns)
  df$RAS <- i$RAS
  df <- gridss.vcftodf.flattenNumeric(df, "RSC", i$RSC, allColumns=allColumns)
  
  df$ASQ <- i$ASQ
  df <- gridss.vcftodf.flattenNumeric(df, "RPQ", i$RPQ, allColumns=allColumns)
  df <- gridss.vcftodf.flattenNumeric(df, "SCQ", i$SCQ, allColumns=allColumns)
  df$RASQ <- i$RASQ
  df <- gridss.vcftodf.flattenNumeric(df, "RSCQ", i$RSCQ, allColumns=allColumns)
  
  df$BAS <- i$BAS
  df <- gridss.vcftodf.flattenNumeric(df, "BRP", i$BRP, allColumns=allColumns)
  df <- gridss.vcftodf.flattenNumeric(df, "BSC", i$BSC, allColumns=allColumns)
  
  df$BASQ <- i$BASQ
  df <- gridss.vcftodf.flattenNumeric(df, "BRPQ", i$BRPQ, allColumns=allColumns)
  df <- gridss.vcftodf.flattenNumeric(df, "BSCQ", i$BSCQ, allColumns=allColumns)
  
  df <- replace(df, is.na(df), 0)
  rownames(df) <- df$variantid
  if (sanityCheck) {
    # should all have mates
    if (!all(df$mate %in% row.names(df))) {
      warning("Unpartnered breakend found")
    }
    mdf <- df[df$mateid,]
    # breakpoint fields should match on both sides of the breakend
    gridss.vcftodf.sanitycheck(df[df$FILTER != mdf$FILTER,], "filter")
    gridss.vcftodf.sanitycheck(df[df$SOMATIC != mdf$SOMATIC,], "somatic")
    gridss.vcftodf.sanitycheck(df[df$IMPRECISE != mdf$IMPRECISE,], "flag")
    gridss.vcftodf.sanitycheck(df[df$SVLEN != mdf$SVLEN,], "svlen")
    gridss.vcftodf.sanitycheck(df[df$SVTYPE != mdf$SVTYPE,], "svtype")
    gridss.vcftodf.sanitycheck(df[df$SPV != mdf$SPV,], "spv")
    gridss.vcftodf.sanitycheck(df[df$QUAL != mdf$QUAL,], "qual")
    gridss.vcftodf.sanitycheck(df[df$CQ != mdf$CQ,], "called qual")
    gridss.vcftodf.sanitycheck(df[df$AS != mdf$RAS,], "assembly")
    gridss.vcftodf.sanitycheck(df[df$RP != mdf$RP,], "read pair")
    gridss.vcftodf.sanitycheck(df[df$SC != mdf$RSC,], "soft clip")
    gridss.vcftodf.sanitycheck(df[abs(df$ASQ - mdf$RASQ) > 0.1,], "assembly qual")
    gridss.vcftodf.sanitycheck(df[abs(df$RPQ - mdf$RPQ) > 0.1,], "read pair qual")
    gridss.vcftodf.sanitycheck(df[abs(df$SCQ - mdf$RSCQ) > 0.1,], "soft clip qual")
  }
  df$hasSC <- paste("SC", ifelse(df$SC > 0 & df$RSC > 0, "Both", ifelse(df$SC > 0, "local", ifelse(df$RSC > 0, "remote", "zero"))))
  df$hasAS <- paste("AS", ifelse(df$AS > 0 & df$RAS > 0, "Both", ifelse(df$AS > 0, "local", ifelse(df$RAS > 0, "remote", "zero"))))
  df$hasRP <- ifelse(df$RP > 0, "Has RP", "No RP")
  df$confidence <- as.factor(ifelse(df$QUAL >= 1000 & df$hasAS == "AS Both", 3, ifelse(df$QUAL >= 500 & df$hasAS != "AS zero", 2, 1)))
  levels(df$confidence) <- c("Low", "Medium", "High")
  df$size <- vcftobpgr(vcf)$size
  return(df)
}
gridss.vcftodf.sanitycheck <- function(failing, desc) {
  failCount <- nrow(failing)
  if (failCount > 0) {
    warning(paste(failCount, "rows failed ", desc, "sanity check "))
  }
}
vcftoroc <- function(vcf, grtp, ...) {
  hits <- breakpointHits(vcftobpgr(vcf), grtp, ...)
  hits$qual <- rowRanges(vcf)$QUAL[hits$queryHits]
  hits <- hits[order(hits$qual),]
  tp <- rep(0, length(grtp))
  tp[hits$subjectHits] <- hits$qual
  tp <- tp[order(-tp)]
  df <- data.frame(qual=tp, sens=1:length(tp) / length(tp))
  df <- df[order(-df$sens),]
  df <- df[!duplicated(df$qual),]
  return(df)
}


