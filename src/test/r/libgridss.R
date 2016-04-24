library(VariantAnnotation)
library(plyr)
source("libvcf.R")

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
  hits$qual <- gridssdf$QUAL[hits[[2]]]
  hits$assemblyboth <- gridssdf$AS[hits[[2]]] > 0 & gridssdf$RAS[hits[[2]]] > 0
  hits$assemblyany <- gridssdf$AS[hits[[2]]] > 0 | gridssdf$RAS[hits[[2]]] > 0
  hits$gridssid <- names(rowRanges(gridssVcf))[hits[[2]]]
  hits$bedid <- names(bed)[hits[[1]]]
  hits <- hits[order(hits$assemblyboth, hits$assemblyany, hits$qual),]
  
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
      df[paste0(name, offset - 1)] <- mat[,offset]
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
  #matchLength <- as.integer(gsub(".*_", "", as.character(i$TRUTH_MATCHES)))
  #mismatchLength <- as.integer(gsub(".*_", "", as.character(i$TRUTH_MISREALIGN)))
  #df$call <- ifelse(!is.na(matchLength), "good", ifelse(!is.na(mismatchLength), "misaligned", "bad"))
  #df$SVLEN <- ifelse(is.na(matchLength), mismatchLength, matchLength)
  df$SVTYPE <- i$SVTYPE
  df$HOMSEQ <- ifelse(is.na(as.character(i$HOMSEQ)), "", as.character(i$HOMSEQ))
  df$HOMLEN <- as.numeric(i$HOMLEN)
  anchorseq <- str_match(as.character(rowRanges(vcf)$ALT), "(.(.*))?[\\[\\]].*[\\[\\]]((.*).)?")[,c(3,5)]
  df$INSSEQ <- ifelse(is.na(anchorseq[,1]), anchorseq[,2], anchorseq[,1])
  df <- gridss.vcftodf.flattenNumeric(df, "REF", i$REF, allColumns=allColumns)
  df <- gridss.vcftodf.flattenNumeric(df, "REFPAIR", i$REFPAIR, allColumns=allColumns)
  df$SPV <- i$SPV
  df$CQ <- i$CQ
  df$BQ <- i$BQ
  
  df$AS <- i$AS
  df <- gridss.vcftodf.flattenNumeric(df, "RP", i$RP, allColumns=allColumns)
  df <- gridss.vcftodf.flattenNumeric(df, "SR", i$SR, allColumns=allColumns)
  df$RAS <- i$RAS
  df <- gridss.vcftodf.flattenNumeric(df, "RSR", i$RSR, allColumns=allColumns)
  df <- gridss.vcftodf.flattenNumeric(df, "ASRP", i$ASRP, allColumns=allColumns)
  df <- gridss.vcftodf.flattenNumeric(df, "ASSR", i$ASSR, allColumns=allColumns)
  df <- gridss.vcftodf.flattenNumeric(df, "ASCRP", i$ASCRP, allColumns=allColumns)
  df <- gridss.vcftodf.flattenNumeric(df, "ASCSR", i$ASCSR, allColumns=allColumns)
  df$ASQ <- i$ASQ
  df <- gridss.vcftodf.flattenNumeric(df, "RPQ", i$RPQ, allColumns=allColumns)
  df <- gridss.vcftodf.flattenNumeric(df, "SRQ", i$SRQ, allColumns=allColumns)
  df$RASQ <- i$RASQ
  df <- gridss.vcftodf.flattenNumeric(df, "RSRQ", i$RSRQ, allColumns=allColumns)
  df$BA <- i$BA
  df <- gridss.vcftodf.flattenNumeric(df, "BUM", i$BUM, allColumns=allColumns)
  df <- gridss.vcftodf.flattenNumeric(df, "BSC", i$BSC, allColumns=allColumns)
  df$BAQ <- i$BAQ
  df <- gridss.vcftodf.flattenNumeric(df, "BUMQ", i$BUMQ, allColumns=allColumns)
  df <- gridss.vcftodf.flattenNumeric(df, "BSCQ", i$BSCQ, allColumns=allColumns)
  
  df <- replace(df, is.na(df), 0)
  rownames(df) <- df$variantid
  if (allColumns) {
    # output aggregate per-category quality scores
    numericSuffix <- 0
    while (length(intersect(c(paste0("RPQ", numericSuffix),
                       paste0("SRQ", numericSuffix),
                       paste0("RSRQ", numericSuffix),
                       paste0("ASCRP", numericSuffix),
                       paste0("ASSR", numericSuffix)),
                       names(df)))) {
      qcolname <- paste0("Q", numericSuffix)
      df[qcolname]  <- 0
      for (prefix in c("RPQ", "SRQ", "RSRQ")) {
        colname <- paste0(prefix, numericSuffix)
        if (colname %in% names(df)) {
          df[qcolname] <- df[qcolname] + df[colname]
        }
      }
      assemblyQual <- (df$ASQ + df$RASQ) * ((df[paste0("ASRP", numericSuffix)] + df[paste0("ASSR", numericSuffix)]) / (df$ASRP + df$ASSR))
      assemblyQual[is.na(assemblyQual)] <- 0
      df[qcolname] <- df[qcolname] + assemblyQual
        
      qcolname <- paste0("Count", numericSuffix)
      df[qcolname]  <- 0
      for (prefix in c("RP", "SR", "RSR")) {
        colname <- paste0(prefix, numericSuffix)
        if (colname %in% names(df)) {
          df[qcolname] <- df[qcolname] + df[colname]
        }
      }
      numericSuffix <- numericSuffix + 1
    }
  }
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
    gridss.vcftodf.sanitycheck(df[df$SR != mdf$RSC,], "soft clip")
    gridss.vcftodf.sanitycheck(df[abs(df$ASQ - mdf$RASQ) > 0.1,], "assembly qual")
    gridss.vcftodf.sanitycheck(df[abs(df$RPQ - mdf$RPQ) > 0.1,], "read pair qual")
    gridss.vcftodf.sanitycheck(df[abs(df$SRQ - mdf$RSCQ) > 0.1,], "soft clip qual")
  }
  df$hasSC <- paste("SC", ifelse(df$SR > 0 & df$RSR > 0, "Both", ifelse(df$SR > 0, "local", ifelse(df$RSR > 0, "remote", "zero"))))
  df$hasAS <- paste("AS", ifelse(df$AS > 0 & df$RAS > 0, "Both", ifelse(df$AS > 0, "local", ifelse(df$RAS > 0, "remote", "zero"))))
  df$hasRP <- ifelse(df$RP > 0, "Has RP", "No RP")
  df$assembly <- ifelse(df$AS > 0 & df$RAS > 0, "Both", ifelse(df$AS + df$RAS > 0, "Single", "No assembly"))
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

