#source("http://bioconductor.org/biocLite.R")
#biocLite("VariantAnnotation", "GenomicFeatures")
#install.packages('reshape')
library(stringr)
library(reshape)
library(parallel)
library(plyr)
library(data.table)
library(VariantAnnotation)
library(GenomicFeatures)

GetMetadataId <- function(filenames) {
  return (regmatches(as.character(filenames), regexpr("[0-9a-f]{32}", as.character(filenames))))
}
GetGenome <- function(file) {
  return (sub(".fa", "", basename(as.character(metadata[GetMetadataId(file),]$CX_REFERENCE))))
}
# Load metadata into a dataframe
LoadMetadata <- function(directory=".") {
  write("Loading metadata", stderr())
  fileList <- list.files(directory, pattern="*.metadata$")
  zeroSizeFiles = file.info(fileList)$size == 0
  if (any(zeroSizeFiles)) {
    warning(paste("Skipping files", fileList[zeroSizeFiles], "as they have 0 size."))
    fileList <- fileList[!zeroSizeFiles]
  }
  #metadata <- foreach (filename=fileList, .export=c("GetMetadataId"), .combine=rbind) %dopar% {
  metadata <- lapply(fileList, function(filename) {
      md <- read.csv(filename, header=FALSE, sep="=", quote = "\"'", col.names=c("CX", "V"))
      md$File <- filename
      md$Id <- GetMetadataId(filename)
      md
  })
  metadata <- do.call(rbind, metadata)
  metadata <- data.frame(lapply(metadata, as.character), stringsAsFactors=FALSE)
  metadata <- cast(metadata, File + Id ~ CX, value="V")  # pivot on context name
  rownames(metadata) <- metadata$Id
  # transform known numeric data to expected type
  metadata$CX_READ_FRAGMENT_LENGTH <- as.numeric(as.character(metadata$CX_READ_FRAGMENT_LENGTH))
  metadata$CX_READ_LENGTH <- as.numeric(as.character(metadata$CX_READ_LENGTH))
  metadata$CX_READ_DEPTH <- as.numeric(as.character(metadata$CX_READ_DEPTH))
  metadata$CX_ALIGNER_SOFTCLIP <- as.numeric(as.character(metadata$CX_ALIGNER_SOFTCLIP))
  write(paste(nrow(metadata), "metadata files loaded"), stderr())
  return(metadata)
}
issymbolic <- function(vcf) {
  if (nrow(vcf) == 0) return(logical(0))
  v <- as.character(unstrsplit(alt(vcf)))
  return(str_detect(v, stringr::fixed("<")) | str_detect(v, stringr::fixed("[")) | str_detect(v, stringr::fixed("]")))
}
svlen <- function(vcf) {
  return(svlentype(vcf)$len)
}
svtype <- function(vcf) {
  return(svlentype(vcf)$type)
}
svlentype <- function(vcf) {
  if (nrow(vcf) == 0) return(list(len=integer(0), type=c()))
  sym <- issymbolic(vcf)
  # calculate from first alt allele
  alleleLen <- width(unlist(alt(vcf)))[c(1, 1 + head(cumsum(elementLengths(alt(vcf))), -1))]
  refLen <- elementLengths(ref(vcf))
  alleleSizeDiff <- alleleLen - refLen
  len=ifelse(sym, NA_integer_, ifelse(alleleLen != refLen, abs(alleleLen - refLen), alleleLen))
  # can't correctly classify events with no length change (INV, DUP, ...) for non-symbolic alleles
  type=ifelse(sym, NA_integer_, ifelse(alleleLen > refLen, "INS", ifelse(alleleLen < refLen, "DEL", "UNKNOWN")))
  
  # override with SV info column data
  svcol <- info(vcf)$SVLEN
  if (!is.null(svcol)) {
    # grab first allele
    len2 <- unlist(svcol)[c(1, 1 + head(cumsum(elementLengths(svcol)), -1))]
    len2[elementLengths(svcol) == 0] <- NA_integer_
    len <- ifelse(is.na(len2), len, len2)
  } else if (!is.null(info(vcf)$END)) {
    # Delly "INFO:SVLEN was a redundant tag because it's just INFO:END - POS for Delly. Since SVLEN is kind of ill-defined in the VCF-Specification I decided to remove it. If you need it for filtering just use end - start."
    len2 <- info(vcf)$END - start(rowRanges(vcf))
    len <- ifelse(is.na(len2), len, len2)
  }
  svcol <- info(vcf)$SVTYPE
  if (!is.null(svcol)) {
    type <- ifelse(is.na(svcol), type, svcol)
  }
  return(list(len=len, type=type))
}
breakendCount <- function(type) {
  return (ifelse(type=="BND", 1,
          ifelse(type %in% c("INS", "DEL", "DUP"), 2,
          ifelse(type=="INV", 4,
          NA_integer_))))
}
# lists overlapping breakpoints
breakpointHits <- function(queryGr, subjectGr, mateQueryGr=queryGr[queryGr$mateIndex,], mateSubjectGr=subjectGr[subjectGr$mateIndex,], ...) {
  if (length(queryGr) != length(mateQueryGr)) {
    stop("queryGr, mateQueryGr have different lengths")
  }
  if (length(queryGr) != length(mateQueryGr)) {
    stop("subjectGr, mateSubjectGr have different lengths")
  }
  dfhits <- rbind(as.data.frame(findOverlaps(queryGr, subjectGr, ...), row.names=NULL),
                  as.data.frame(findOverlaps(mateQueryGr, mateSubjectGr, ...), row.names=NULL))
  dfhits <- dfhits[duplicated(dfhits),] # both breakends match
  return(dfhits)
}
# counts overlapping breakpoints
countBreakpointHits <- function(queryGr, subjectGr, mateQueryGr=queryGr[queryGr$mate,], mateSubjectGr=subjectGr[subjectGr$mate,], ...) {
  dfhits <- breakpointHits(queryGr, subjectGr, mateQueryGr, mateSubjectGr, ...)
  queryHits <- rep(0, length(queryGr))
  queryHits[count(dfhits, "queryHits")$queryHits] <- count(dfhits, "queryHits")$freq
  subjectHits <- rep(0, length(subjectGr))
  subjectHits[count(dfhits, "subjectHits")$subjectHits] <- count(dfhits, "subjectHits")$freq
  return(list(queryHitCount=queryHits, subjectHitCount=subjectHits))
}
# converts a VCF to a GRanges containing the paired variant breakends with the following fields:
# vcfIndex: index of variant in VCF
# mateIndex: index of matching breakend in the GRanges object
# size: event size
# SVTYPE: type of event called
vcftobpgr <- function(vcf) {
  lentype <- svlentype(vcf)
  grcall <- rowRanges(vcf)
  grcall$vcfIndex <- seq_along(grcall)
  grcall$mateIndex <- rep(NA_integer_, length(grcall))
  grcall$SVTYPE <- lentype$type
  grcall$size <- lentype$len
  grcall$untemplated <- rep(NA_integer_, length(grcall))
  strand(grcall) <- rep("+", length(grcall))
  rows <- !issymbolic(vcf)
  if (!is.null(grcall$SVTYPE)) {
    rows <- rows & is.na(grcall$SVTYPE)
  }
  grcall$cistartoffset <- rep(0, length(grcall))
  grcall$ciwidth <- rep(0, length(grcall))
  if ("CIPOS" %in% names(info(vcf))) {
    # Expand call position by CIPOS
    offsets <- matrix(unlist(info(vcf)$CIPOS), ncol = 2, byrow = TRUE)
    offsets[is.na(offsets)] <- 0
    grcall$cistartoffset <- offsets[,1]
    grcall$ciwidth <- offsets[,2] - offsets[,1]
  }
  if (any(rows)) {
    # non-symbolic VCF record
    browser()
    stop("TODO: handle non-symbolic alleles")
  }
  rows <- grcall$SVTYPE=="BND"
  if (any(rows)) {
    # set strand for BND
    bndMatches <- str_match(as.character(rowRanges(vcf[rows,])$ALT), "(.*)(\\[|])(.*)(\\[|])(.*)")
    preBases <- bndMatches[,2]
    bracket <- bndMatches[,3]
    remoteLocation <- bndMatches[,4]
    postBases <- bndMatches[,6]
    strand(grcall[rows,]) <- ifelse(str_length(preBases) > 0, "+", "-")
    if (!is.null(info(vcf)$MATEID)) {
      grcall[rows,]$mateIndex <- match(as.character(info(vcf)$MATEID[rows]), names(rowRanges(vcf)))
    } else if (!is.null(info(vcf)$PARID)) {
      grcall[rows,]$mateIndex <- match(info(vcf)$PARID[rows], names(rowRanges(vcf)))
    }
    grcall[rows,]$untemplated <- str_length(preBases) + str_length(postBases) - str_length(as.character(rowRanges(vcf)$REF)[rows])
    
    mateBnd <- grcall[grcall[rows,]$mateIndex,]
    grcall[rows,]$size <- ifelse(seqnames(mateBnd)==seqnames(grcall[rows,]), abs(start(grcall[rows,]) - start(mateBnd)) + grcall[rows,]$untemplated, NA_integer_)
  }
  # non-standard event type used by DELLY
  rows <- grcall$SVTYPE=="TRA"
  if (any(rows)) {
    strand(grcall[rows,]) <- ifelse(substr(info(vcf)$CT[rows], 1, 1) == "3", "+", "-")
    grcall[rows]$mateIndex <- length(grcall) + seq_len(sum(rows))
    eventgr <- grcall[rows]
    eventgr$mateIndex <- seq_len(length(grcall))[rows]
    # need to make new copy since seqlevels might not match
    bareeventgr <- GRanges(
      seqnames=info(vcf)$CHR2[rows],
      ranges=IRanges(start=info(vcf)$END[rows], width=1),
      strand=ifelse(substr(info(vcf)$CT[rows], 4, 4) == "3", "+", "-"))
    mcols(bareeventgr) <- mcols(eventgr)
    grcall <- c(grcall, eventgr)
  }
  rows <- grcall$SVTYPE=="DEL"
  if (any(rows)) {
    grcall[rows]$mateIndex <- length(grcall) + seq_len(sum(rows))
    eventgr <- grcall[rows]
    strand(eventgr) <- "-"
    ranges(eventgr) <- IRanges(start=start(eventgr) + abs(grcall$size[rows]), width=1)
    eventgr$mateIndex <- seq_len(length(grcall))[rows]
    grcall <- c(grcall, eventgr)
  }
  rows <- grcall$SVTYPE=="INS"
  if (any(rows)) {
    grcall[rows]$mateIndex <- length(grcall) + seq_len(sum(rows))
    eventgr <- grcall[rows]
    strand(eventgr) <- "-"
    ranges(eventgr) <-IRanges(start=start(eventgr) + 1, width=1)
    eventgr$mateIndex <- seq_len(length(grcall))[rows]
    grcall <- c(grcall, eventgr)
  }
  rows <- grcall$SVTYPE=="INV"
  if (any(rows)) {
    grcall[rows]$mateIndex <- length(grcall) + seq_len(sum(rows))
    eventgr1 <- grcall[rows]
    eventgr2 <- grcall[rows]
    eventgr3 <- grcall[rows]
    strand(eventgr2) <- "-"
    strand(eventgr3) <- "-"
    ranges(eventgr1) <- IRanges(start=start(eventgr1) + abs(grcall[rows]$size), width=1)
    ranges(eventgr3) <- IRanges(start=start(eventgr3) + abs(grcall[rows]$size), width=1)
    eventgr1$mateIndex <- seq_len(length(grcall))[rows]
    eventgr2$mateIndex <- length(grcall) + length(eventgr1) + length(eventgr2) + seq_len(length(eventgr3))
    eventgr3$mateIndex <- length(grcall) + length(eventgr1) + seq_len(length(eventgr2))
    grcall <- c(grcall, eventgr1, eventgr2, eventgr3)
  }
  rows <- grcall$SVTYPE=="DUP" | grcall$SVTYPE=="DUP:TANDEM" 
  if (any(rows)) {
    # note: pindel tandem duplication includes the sequence being duplicated in the REF allele
    grcall[rows]$mateIndex <- length(grcall) + seq_len(sum(rows))
    eventgr <- grcall[rows]
    strand(grcall[rows]) <- "-"
    ranges(eventgr) <- IRanges(start=start(eventgr) + abs(grcall[rows]$size), width=1)
    eventgr$mateIndex <- seq_len(length(grcall))[rows]
    grcall <- c(grcall, eventgr)
  }
  rows <- grcall$SVTYPE=="RPL"
  if (any(rows)) {
    # pindel 'replacement' SV type for handling untemplated sequence
    # place breakpoints at start and end of ref allele position
    # (not ideal, but it's the best conversion we can do without examine the actual sequences)
    grcall[rows]$mateIndex <- length(grcall) + seq_len(sum(rows))
    eventgr <- grcall[rows]
    strand(eventgr) <- "-"
    ranges(eventgr) <- IRanges(start=start(eventgr) + abs(elementLengths(ref(vcf))[rows]), width=1)
    eventgr$mateIndex <- seq_len(length(grcall))[rows]
    grcall <- c(grcall, eventgr)
  }
  width(grcall) <- rep(1, length(grcall))
  start(grcall) <- start(grcall) + grcall$cistartoffset
  end(grcall) <- start(grcall) + grcall$ciwidth
  if (any(is.na(grcall$mateIndex))) {
    browser()
    stop(paste0("Unhandled SVTYPE ", unique(grcall$SVTYPE[is.na(grcall$mateIndex)])))
  }
  # Check the partner of the partner of each row is the row itself
  if (!all(grcall[grcall$mateIndex,]$mateIndex == seq_along(grcall))) {
    browser()
    stop("Breakends are not uniquely paired.")
  }
  if (any(grcall$mateIndex == seq_along(grcall))) {
    browser()
    stop("Breakend cannot be partner of itself - have all appropriate SV types been handled?")
  }
  return(grcall)
}
interval_distance <- function(s1, e1, s2, e2) {
  return (ifelse(s2 >= s1 & s2 <= e1, 0,
          ifelse(s1 >= s2 & s1 <= e2, 0,
          ifelse(s1 < s2, s2 - e1, s1 - e2))))
}
distanceToClosest <- function(query, subject) {
  distanceHits <- distanceToNearest(query, subject)
  result <- rep(NA, length(query))
  result[queryHits(distanceHits)] <- as.data.frame(distanceHits)$distance
  return(result)
}


CalculateTruth <- function(callvcf, truthvcf, maxerrorbp, maxerrorpercent=0.25, ...) {
  if (any(!is.na(rowRanges(callvcf)$QUAL) & rowRanges(callvcf)$QUAL < 0)) {
    stop("Precondition failure: variant exists with negative quality score")
  }
  grcall <- vcftobpgr(callvcf)
  grtruth <- vcftobpgr(truthvcf)
  #TMPDEBUG: grtruth[overlapsAny(grtruth, GRanges("chr12", IRanges(148000, width=1000))),]
  hits <- breakpointHits(query=grcall, subject=grtruth, maxgap=maxerrorbp + max(0, grcall$ciwidth), ...)
  hits$QUAL <- grcall$QUAL[hits$queryHits]
  hits$QUAL[is.na(hits$QUAL)] <- 0
  hits$poserror <- interval_distance(start(grcall)[hits$queryHits], end(grcall)[hits$queryHits], start(grtruth)[hits$subjectHits], end(grtruth)[hits$subjectHits]) # breakend error distribution
  hits$calledsize <- abs(grcall$size[hits$queryHits])
  hits$expectedsize <- abs(grtruth$size[hits$subjectHits])
  hits$errorsize <- abs(abs(hits$calledsize - hits$expectedsize) - grcall$ciwidth[hits$queryHits])
  hits$percentsize <- hits$calledsize / hits$expectedsize
  # TODO: add untemplated into sizerror calculation for insertions
  hits <- hits[is.na(hits$poserror) | hits$poserror <= maxerrorbp, ] # position must be within margin
  hits <- hits[is.na(hits$expectedsize) | is.na(hits$calledsize) | (hits$percentsize >= 1 - maxerrorpercent & hits$percentsize <= 1 + maxerrorpercent), ] # size must approximately match
  # TODO: filter mismatched event types (eg: DEL called for INS)
  
  # per truth variant
  hits <- data.table(hits, key="subjectHits")
  tdf <- hits[, list(hits=.N, poserror=min(poserror), errorsize=min(errorsize), QUAL=max(QUAL)), by="subjectHits"] 
  grtruth$QUAL <- -1
  grtruth$QUAL[tdf$subjectHits] <- tdf$QUAL
  grtruth$hits <- 0
  grtruth$hits[tdf$subjectHits] <- tdf$hits
  grtruth$poserror <- NA_integer_
  grtruth$poserror[tdf$subjectHits] <- tdf$poserror
  grtruth$errorsize <- NA_integer_
  grtruth$errorsize[tdf$subjectHits] <- tdf$errorsize
  grtruth$tp <- grtruth$hits > 0
  truthdf <- data.frame(
    vcfIndex=seq_along(rowRanges(truthvcf)),
    SVTYPE=svtype(truthvcf),
    SVLEN=svlen(truthvcf),
    expectedbehits=breakendCount(svtype(truthvcf))
    )
  tdfbe <- data.table(as.data.frame(mcols(grtruth)), key="vcfIndex")[, list(
    behits=sum(tp),
    poserror=mean(poserror),
    maxposerror=max(poserror),
    errorsize=min(errorsize),
    QUAL=mean(QUAL)
    ), by="vcfIndex"]
  if (any(truthdf$vcfIndex != tdfbe$vcfIndex)) {
    stop("Sanity check failure: tdfbe vcfIndex offsets do not match truthdf")
  }
  tdfbe$vcfIndex <- NULL
  truthdf <- cbind(truthdf, tdfbe)
  truthdf$tp <- truthdf$expectedbehits == truthdf$behits
  truthdf$partialtp <- !truthdf$tp & truthdf$behits > 0
  truthdf$fn <- !truthdf$tp
  
  # per variant call
  hits <- data.table(hits, key="queryHits")
  cdf <- hits[, list(hits=.N, poserror=min(poserror), errorsize=min(errorsize)), by="queryHits"]
  grcall$hits <- rep(0, length(grcall))
  grcall$hits[cdf$queryHits] <- cdf$hits
  grcall$poserror <- rep(NA_integer_, length(grcall))
  grcall$poserror[cdf$queryHits] <- cdf$poserror
  grcall$errorsize <- rep(NA_integer_, length(grcall))
  grcall$errorsize[cdf$queryHits] <- cdf$errorsize
  grcall$tp <- grcall$hits > 0
  calldf <- data.frame(
    vcfIndex=seq_along(rowRanges(callvcf)),
    QUAL=rowRanges(callvcf)$QUAL,
    SVTYPE=svtype(callvcf),
    SVLEN=svlen(callvcf),
    expectedbehits=breakendCount(svtype(callvcf))
  )
  cdfbe <- data.table(as.data.frame(mcols(grcall)), key="vcfIndex")[, list(
    behits=sum(tp),
    poserror=mean(poserror),
    errorsize=min(errorsize),
    maxposerror=max(poserror)
    ), by="vcfIndex"]
  if (any(calldf$vcfIndex != cdfbe$vcfIndex)) {
    stop("Sanity check failure: cdfbe vcfIndex offsets do not match calldf")
  }
  cdfbe$vcfIndex <- NULL
  calldf <- cbind(calldf, cdfbe)
  calldf$tp <- calldf$expectedbehits == calldf$behits
  calldf$partialtp <- !calldf$tp & calldf$behits > 0
  calldf$fp <- !calldf$tp
  return(list(calls=calldf, truth=truthdf))
}
CalculateTruthSummary <- function(vcfs, maxerrorbp, ...) {
  truthSet <- lapply(vcfs, function(vcf) {
    md <- attr(vcf, "metadata")
    if (is.null(md)) {
      warning("VCF missing metadata - ignoring")
      return(NULL)
    }
    if (is.null(md$CX_CALLER) || is.na(md$CX_CALLER)) {
      return(NULL)
    }
    if (is.null(md$CX_REFERENCE_VCF) || is.na(md$CX_REFERENCE_VCF)) {
      warning(paste0(md$Id, " missing reference vcf"))
      return(NULL)
    }
    truthvcf <- vcfs[[GetMetadataId(md$CX_REFERENCE_VCF)]]
    if (is.null(truthvcf)) {
      warning(paste0("Missing truth for", md$Id))
      return(NULL)
    }
    write(paste0("Loading truth for ", md$Id), stderr())
    callTruthPair <- CalculateTruth(vcf, truthvcf, maxerrorbp, ...)
    if (nrow(callTruthPair$calls) > 0) {
      callTruthPair$calls <- cbind(callTruthPair$calls, md, row.names=NULL)
    }
    callTruthPair$truth <- cbind(callTruthPair$truth, md, row.names=NULL)
    return(callTruthPair)
  })
  truthSet[sapply(truthSet, is.null)] <- NULL # Remove NULLs
  calls = rbindlist(lapply(truthSet, function(x) x$calls), use.names=TRUE, fill=TRUE)
  truth = rbindlist(lapply(truthSet, function(x) x$truth), use.names=TRUE, fill=TRUE)
  return(list(calls=calls, truth=truth))
}
TruthSummaryToROC <- function(ts, bylist=c("CX_ALIGNER", "CX_ALIGNER_SOFTCLIP", "CX_CALLER", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_LENGTH", "CX_REFERENCE_VCF_VARIANTS", "SVTYPE")) {
  callset <- ts$calls[, c(bylist, "QUAL", "tp", "fp"), with=FALSE]
  callset$fn <- FALSE
  callset$QUAL[is.na(callset$QUAL)] <- 0
  callset[callset$tp==FALSE,] # TPs are in both call and truth sets - drop the call set version
  truthset <- ts$truth[, c(bylist, "QUAL", "tp", "fn"), with=FALSE]
  truthset$fp <- FALSE
  truthset$QUAL[is.na(truthset$QUAL)] <- -1
  combined <- rbind(callset, truthset)
  combined$tp <- as.integer(combined$tp)
  combined$fn <- as.integer(combined$fn)
  combined$fp <- as.integer(combined$fp)
  setkeyv(combined, bylist)
  # aggregate from high to low
  combined <- combined[order(-combined$QUAL),]
  combined[,`:=`(ntruth=sum(tp)+sum(fn), ncalls=sum(tp)+sum(fp), tp=cumsum(tp), fp=cumsum(fp), fn=cumsum(fn), QUAL=cummin(QUAL)), by=bylist]
  combined <- combined[!duplicated(combined[, c(bylist, "QUAL"), with=FALSE], fromLast=TRUE),] # take only one data point per QUAL
  combined[,`:=`(sens=tp/ntruth, prec=tp/(tp+fp), fdr=fp/(tp+fp))]
  return(combined)
}
FilterOutSNV <- function(vcf, caller) {
  if (any(elementLengths(alt(z)) != 1)) {
    browser()
    stop(paste("Analysis not designed for multiple alleles in a single VCF record.", filename))
  }
  # TODO: remove SNPs
  # - same length
  # - differ by 1 base (samtools will sometimes call a SNP as AAAC vs AAAT)
  # - not symbolic (no < or [ )
  #if (caller <- c("gatk")) {
  #  # strip SNPs based on alt allele length matching ref allele
  #  vcf <- vcf[elementLengths(stringSet)]
  #  vcf <- vcf[unlist(lapply(alt(vcf), function (stringSet) { elementLengths(stringSet)[1] } )) != elementLengths(ref(vcf))]
  #}
  return(vcf)
}
# Load VCFs into a list
LoadVcfs <- function(metadata, directory=".", pattern="*.vcf$") {
  write("Loading VCFs", stderr())
  fileList <- list.files(directory, pattern=pattern)
  zeroSizeFiles = file.info(fileList)$size == 0
  if (any(zeroSizeFiles)) {
    write(paste("Skipping file", fileList[zeroSizeFiles], "due to 0 size.\n"))
    warning(paste("Skipping files", paste(fileList[zeroSizeFiles]), "due to 0 size.\n"))
    fileList <- fileList[!zeroSizeFiles]
  }
  # Parallel load of VCFs
  #vcfs <- foreach (filename=fileList, .packages="VariantAnnotation") %dopar% {
  vcfs <- lapply(fileList, function(filename) {
    write(paste0("Loading ", filename), stderr())
    vcf <- readVcf(filename, "unknown")
    attr(vcf, "filename") <- filename
    vcf
  })
  vcfs <- lapply(vcfs, function(vcf) {
    filename <- attr(vcf, "filename")
    id <- GetMetadataId(filename)
    md <- metadata[id, ]
    attr(vcf, "id") <- id
    attr(vcf, "metadata") <- md
    vcf
  })
  names(vcfs) <- GetMetadataId(fileList)
  vcfs[sapply(vcfs, is.null)] <- NULL # Remove NULL VCFs list
  write(paste("Loaded", length(vcfs), "VCFs"), stderr())
  return(vcfs)
}
