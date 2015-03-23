library(VariantAnnotation)
library(plyr)


gridss.annotateBreakpointHits <- function(bed, bedMate, gridssVcf, ...) {
  # filter breakends and one-sided breakpoint calls
  gridssVcf <- gridssVcf[as.character(info(gridssVcf)$MATEID) %in% row.names(gridssVcf),]
  gridssdf <- gridss.vcftodf(gridssVcf)
  callPos <- rowData(gridssVcf)
  callPos$mate <- as.character(info(gridssVcf)$MATEID)
  strand(callPos) <- ifelse(str_detect(as.character(callPos$ALT), "[[:alpha:]]+(\\[|]).*(\\[|])"), "+", "-")
  callPosMate <- callPos[callPos$mate,]
  
  hits <- rbind(
    as.data.frame(findOverlaps(bed, callPos, ...)),
    as.data.frame(findOverlaps(bed, callPosMate, ...)),
    as.data.frame(findOverlaps(bedMate, callPos, ...)),
    as.data.frame(findOverlaps(bedMate, callPosMate, ...)))
  hits <- hits[duplicated(hits),]
  hits$qual <- rowData(gridssVcf)$QUAL[hits[[2]]]
  hits$hc <- gridssdf$QUAL[hits[[2]]] >= 1000 & gridssdf$AS[hits[[2]]] > 0 & gridssdf$RAS[hits[[2]]] > 0
  annotatedBed <- bed
  annotatedBed$called <- "miss"
  annotatedBed$called[hits[[1]]] <- "hit"
  annotatedBed$called[hits[hits$hc,][[1]]] <- "high confidence"
  annotatedBed$qual <- 0
  annotatedBed$qual[aggregate(hits, FUN=max, by=list(hits$queryHits))$queryHits] <- aggregate(hits, FUN=max, by=list(hits$queryHits))$qual
  return(annotatedBed)
}

gridss.overlaps <- function(vcf, bed, ...) {
  hits <- findOverlaps(rowData(vcf), bed, ...)
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
# adds tumour/normal/total columns to the given data frame
gridss.truthdetails.processvcf.addtn <- function(df, name, column, max=FALSE) {
  if (max) {
    df[[paste0(name)]] <- pmax(ciltoarray(column, 1), ciltoarray(column, 2))
  } else {
    df[[paste0(name)]] <- ciltoarray(column, 1) + ciltoarray(column, 2)
  }
  df[[paste0(name, "Normal")]] <- ciltoarray(column, 1)
  df[[paste0(name, "Tumour")]] <- ciltoarray(column, 2)
  
  is.na(df[[paste0(name)]]) <- 0
  is.na(df[[paste0(name, "Normal")]]) <- 0
  is.na(df[[paste0(name, "Tumour")]]) <- 0
  return(df)
}
gridss.removeUnpartnerededBreakend <- function(vcf) {
  toKeep <- as.character(info(vcf)$MATEID) %in% row.names(vcf)
  if (sum(toKeep == FALSE) > 0) {
    warning(paste("Found ", sum(toKeep == FALSE), " unpartnered breakends."))
  }
  return(vcf[toKeep,])
}
gridss.vcftodf <- function(vcf, sanityCheck=TRUE) {
  i <- info(vcf)
  df <- data.frame(variantid=names(rowData(vcf)))
  df$FILTER=rowData(vcf)$FILTER
  df$QUAL <- rowData(vcf)$QUAL
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
  df <- gridss.truthdetails.processvcf.addtn(df, "REF", i$REF)
  df <- gridss.truthdetails.processvcf.addtn(df, "REFPAIR", i$REFPAIR)
  df$SPV <- i$SPV
  df$CQ <- i$CQ
  df$BQ <- i$BQ
  
  df$AS <- i$AS
  df <- gridss.truthdetails.processvcf.addtn(df, "RP", i$RP)
  df <- gridss.truthdetails.processvcf.addtn(df, "SC", i$SC)
  df$RAS <- i$RAS
  df <- gridss.truthdetails.processvcf.addtn(df, "RSC", i$RSC)
  
  df$ASQ <- i$ASQ
  df <- gridss.truthdetails.processvcf.addtn(df, "RPQ", i$RPQ)
  df <- gridss.truthdetails.processvcf.addtn(df, "SCQ", i$SCQ)
  df$RASQ <- i$RASQ
  df <- gridss.truthdetails.processvcf.addtn(df, "RSCQ", i$RSCQ)
  
  df$BAS <- i$BAS
  df <- gridss.truthdetails.processvcf.addtn(df, "BRP", i$BRP)
  df <- gridss.truthdetails.processvcf.addtn(df, "BSC", i$BSC)
  
  df$BASQ <- i$BASQ
  df <- gridss.truthdetails.processvcf.addtn(df, "BRPQ", i$BRPQ)
  df <- gridss.truthdetails.processvcf.addtn(df, "BSCQ", i$BSCQ)
  
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
  hits$qual <- rowData(vcf)$QUAL[hits$queryHits]
  hits <- hits[order(hits$qual),]
  tp <- rep(0, length(grtp))
  tp[hits$subjectHits] <- hits$qual
  tp <- tp[order(-tp)]
  df <- data.frame(qual=tp, sens=1:length(tp) / length(tp))
  df <- df[order(-df$sens),]
  df <- df[!duplicated(df$qual),]
  return(df)
}


