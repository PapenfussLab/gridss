#source("http://bioconductor.org/biocLite.R")
#biocLite("VariantAnnotation", "GenomicFeatures")
#install.packages('reshape')

library(VariantAnnotation)
library(GenomicFeatures)
library(reshape)
library(parallel)
library(plyr)

GetMetadataId <- function(filenames) {
  return (regmatches(as.character(filenames), regexpr("[0-9a-f]{32}", as.character(filenames))))
}
GetGenome <- function(file) {
  return (sub(".fa", "", basename(as.character(metadata[GetMetadataId(file),]$CX_REFERENCE))))
}
# Load metadata into a dataframe
LoadMetadata <- function() {
  write("Loading metadata", stderr())
  metadata <- lapply(list.files(".", pattern="*.metadata$"), function(filename) {
    tmp <- read.csv(filename, header=FALSE, sep="=", quote = "\"'", col.names=c("CX", "V"))
    tmp$File <- filename
    tmp$Id <- GetMetadataId(filename)
    return (tmp)
  })
  metadata <- do.call("rbind", metadata)  # collapse data frames # stringsAsFactors
  metadata <- data.frame(rapply(metadata, as.character, classes="factor", how="replace"), stringsAsFactors=FALSE)  # workaround for do.call("rbind" not liking stringsAsFactors=FALSE
  #metadata <- melt(metadata)
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
svlen <- function(vcf) {
  if (nrow(vcf) == 0) return(integer(0))
  svcol <- info(vcf)$SVLEN
  if (!is.null(svcol)) {
    # grab first allele
    len <- unlist(svcol)[c(1, 1 + head(cumsum(elementLengths(svcol)), -1))]
    len[elementLengths(svcol) == 0] <- NA
    return(len)
  }
  if (attr(vcf, "metadata")$CX_CALLER %in% c("samtools", "gatk", "dindel", "SOAPindel", "soapindel")) {
    # calculate from genotype
    # grab the length of the first alt allele
    # Notes:
    # unlist flattens multiple alleles
    # cumsum of elementLenths returns offset of the last allele in the flattened list
    # to get the first, we shift over by 1 so [1, 3, 4] becomes [+1+, 2, 4, -5-]
    return(width(unlist(alt(vcf))[c(1, 1 + head(cumsum(elementLengths(alt(vcf))), -1))])  - elementLengths(ref(vcf)))
    #return(unlist(lapply(alt(vcf), function (stringSet) { elementLengths(stringSet)[1] } )) - elementLengths(ref(vcf)))
  }
  stop("Unable to determine SVLEN from INFO columns")
}
svtype <- function(vcf) {
  svcol <- info(vcf)$SVTYPE
  if (!is.null(svcol)) {
    return(svcol)
  } else if (length(fixed(vcf)) == 0) {
    return(c())
  } else {
    # fall-back to infer from SVLEN
    return(ifelse(VcfGet$SVLEN(vcf) > 0, "INS", "DEL"))
  }
}
# lists overlapping breakpoints
breakpointHits <- function(queryGr, subjectGr, mateQueryGr=queryGr[queryGr$mate,], mateSubjectGr=subjectGr[subjectGr$mate,], ...) {
  dfhits <- rbind(as.data.frame(findOverlaps(queryGr, subjectGr, ...)),
                  as.data.frame(findOverlaps(mateQueryGr, mateSubjectGr, ...)))
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
# counts overlapping breakpoints
# VCF records MUST each be a valid symbolic breakend alleles
# VCF MUST have a single valid MATEID for every record
# GRanges MUST have a $mate set to the name of the other breakend for that breakpoint
countVcfGrBreakpointHits <- function(vcf, gr, ...) {
  vcfgr <- vcftobpgr(vcf)
  return(countBreakpointHits(vcfgr, gr, vcfgr[vcfgr$mate,], gr[gr$mate,], ...))
}
# converts a VCF to a GRanges containing the paired variant breakends with the following fields:
# vcfIndex: index of variant in VCF
# mateIndex: index of matching breakend in the GRanges object
# size: event size
# SVTYPE: type of event called
vcftobpgr <- function(vcf) {
  vcfgr <- rowData(vcf)
  vcfgr$vcfIndex <- seq(1:nrow(vcf))
  vcfgr$mateIndex <- NA_integer_ # seq(1:nrow(vcf))
  vcfgr$SVTYPE <- info(vcf)$SVTYPE
  vcfgr$size <- 0
  vcfgr$size <- svlen(vcf)
  strand(vcfgr) <- "+"
  if (any(vcfgr$SVTYPE=="BND")) {
    # set strand for BND
    strand(vcfgr[vcfgr$SVTYPE=="BND",]) <- ifelse(str_detect(as.character(rowData(vcf[vcfgr$SVTYPE=="BND",])$ALT), "[[:alpha:]]+(\\[|]).*(\\[|])"), "+", "-")
    vcfgr$mate <- NULL
    if (!is.null(info(vcf)$MATEID)) {
      vcfgr[vcfgr$SVTYPE=="BND",]$mateIndex <- match(as.character(info(vcf)$MATEID[vcfgr$SVTYPE=="BND"]), names(rowData(vcf)))
    } else if (!is.null(info(vcf)$PARID)) {
      vcfgr[vcfgr$SVTYPE=="BND",]$mateIndex <- match(info(vcf)$PARID[vcfgr$SVTYPE=="BND"], names(rowData(vcf)))
    }
  }
  if (any(vcfgr$SVTYPE=="DEL")) {
    vcfgr[vcfgr$SVTYPE=="DEL"]$mateIndex <- length(vcfgr) + seq(1, sum(vcfgr$SVTYPE=="DEL"))
    delgr <- vcfgr[vcfgr$SVTYPE=="DEL"]
    strand(delgr) <- "-"
    start(delgr) <- start(delgr) - vcfgr[vcfgr$SVTYPE=="DEL"]$size
    delgr$mateIndex <- seq(1, length(vcfgr))[vcfgr$SVTYPE=="DEL"]
    vcfgr <- c(vcfgr, delgr)
  }
  if (any(vcfgr$SVTYPE=="INS")) {
    vcfgr[vcfgr$SVTYPE=="INS"]$mateIndex <- length(vcfgr) + seq(1, sum(vcfgr$SVTYPE=="INS"))
    insgr <- vcfgr[vcfgr$SVTYPE=="INS"]
    strand(insgr) <- "-"
    ranges(insgr) <-IRanges(start=start(insgr) + 1, width=1)
    insgr$mateIndex <- seq(1, length(vcfgr))[vcfgr$SVTYPE=="INS"]
    vcfgr <- c(vcfgr, insgr)
  }
  if (any(vcfgr$SVTYPE=="INV")) {
    vcfgr[vcfgr$SVTYPE=="INV"]$mateIndex <- length(vcfgr) + seq(1, sum(vcfgr$SVTYPE=="INV"))
    inv1gr <- vcfgr[vcfgr$SVTYPE=="INV"]
    inv2gr <- vcfgr[vcfgr$SVTYPE=="INV"]
    inv3gr <- vcfgr[vcfgr$SVTYPE=="INV"]
    strand(inv2gr) <- "-"
    strand(inv3gr) <- "-"
    ranges(inv1gr) <-IRanges(start=start(inv1gr) + abs(vcfgr[vcfgr$SVTYPE=="INV"]$size), width=1)
    ranges(inv3gr) <-IRanges(start=start(inv3gr) + abs(vcfgr[vcfgr$SVTYPE=="INV"]$size), width=1)
    inv1gr$mateIndex <- seq(1, length(vcfgr))[vcfgr$SVTYPE=="INV"]
    inv2gr$mateIndex <- length(vcfgr) + length(inv1gr) + length(inv2gr) + seq(1, length(inv3gr))
    inv3gr$mateIndex <- length(vcfgr) + length(inv1gr) + seq(1, length(inv2gr))
    vcfgr <- c(vcfgr, inv1gr, inv2gr, inv3gr)
  }
  width(ranges(vcfgr)) <- 1
  # TODO handle CIPOS by widening the breakend width
  if (FALSE & "CIPOS" %in% names(info(vcf))) {
    # Expand call position by CIPOS
    offsets <- matrix(unlist(info(vcf)$CIPOS), ncol = 2, byrow = TRUE)
    offsets[is.na(offsets)] <- 0
    start(ranges(vcfgr)) <- start(ranges(vcfgr)) + offsets[,1]
    end(ranges(vcfgr)) <- end(ranges(vcfgr)) + offsets[,2]
  }
  # Check the partner of the partner of each row is the row itself
  if (!all(vcfgr[vcfgr$mateIndex,]$mateIndex == seq(1:length(vcfgr)))) {
    stop("Breakends are not uniquely paired.")
  }
  if (any(vcfgr$mateIndex == seq(1:length(vcfgr)))) {
    stop("Breakend cannot be partner of itself - have all appropriate SV types been handled?")
  }
  return(vcfgr)
}

# Load VCFs into a list
LoadVcfs <- function(metadata, pattern="[0-9a-f]{32}(.reference)?.vcf$") {
  write("Loading VCFs", stderr())
  vcfFileList <- list.files(".", pattern=pattern)
  LoadVcf <- function(filename) {
    write(paste("Loading", filename), stderr())
    id <- GetMetadataId(filename)
    md <- metadata[id, ]
    if (md$CX_CALLER %in% c("socrates", "hydra")) {
      # for now, skip socrates and hydra since readVcf fails with only breakend-only VCFs
      return (NULL)
    }
    vcf <- readVcf(filename, GetGenome(filename))
    if (is.null(vcf)) error(paste("Unable to parse", filename))
    attr(vcf, "id") <- id
    attr(vcf, "metadata") <- md
    if (nrow(vcf) == 0) {
      warning(paste(filename, "has 0 records. Check", md$CX_CALLER, "script is working."))
    } else if (!is.null(info(vcf)$SVTYPE)) {
      # strip SNPs
      vcf <- vcf[which(!is.na(info(vcf)$SVTYPE))]
      if (nrow(vcf) == 0) {
        warning(paste(filename, "has 0 SV records. Check SVTYPE field is populated."))
      }
    } else if (!is.null(info(vcf)$INDEL)) {
      # strip SNPs from samtools output
      vcf <- vcf[info(vcf)$INDEL]
      if (length(fixed(vcf)) == 0) {
        warning(paste(filename, "has 0 INDEL records."))
      }
    }
    if (any(sapply(fixed(vcf)$ALT, length) != 1)) {
      if (md$CX_CALLER %in% c("samtools", "gatk", "dindel")) {
        warning("Calling only first alternate allele for each SV record")
      } else {
        browser()
        stop(paste("Analysis not designed for multiple alleles in a single INDEL/SV VCF record.", filename))
      }
    }
    if (md$CX_CALLER %in% c("gatk")) {
      # strip SNPs based on alt allele length matching ref allele
      vcf <- vcf[unlist(lapply(alt(vcf), function (stringSet) { elementLengths(stringSet)[1] } )) != elementLengths(ref(vcf))]
    }
    if (nrow(vcf) > 0 && (md$CX_CALLER %in% c("clever", "breakdancer-max", "gasvpro", "delly", "crest") || all(rowData(vcf)@ranges@width == 1))) {
      # fix broken IRanges (0 length, or all 1) by using SVLEN field
      if (!is.null(info(vcf)$SVLEN)) {
        # delete = 1 + size of deletion
        # insert = 1
        if (!is.null(info(vcf)$SVTYPE)) {
          rowData(vcf)@ranges@width <- ifelse(
            # don't change BND
            info(vcf)$SVTYPE=="BND",
            rowData(vcf)@ranges@width,
            as.integer(rowData(vcf)@ranges@width + ifelse(as.integer(info(vcf)$SVLEN) > 0, 1, 1 - as.integer(info(vcf)$SVLEN))))  
        } else {
          rowData(vcf)@ranges@width <- as.integer(rowData(vcf)@ranges@width + ifelse(as.integer(info(vcf)$SVLEN) > 0, 1, 1 - as.integer(info(vcf)$SVLEN)))
        }
      }
    }
    if (any(is.na(rowData(vcf)@ranges@width))) {
      browser()
      stop("Sanity check failure: NA width")
    }
    # for now, skip indel sanity checks
    return (vcf)
    # Sanity checks
    if (any(rowData(vcf)@ranges@width <= 0)) {
      browser()
      # findOverlaps asymetrically breaks with 0 width IRanges
      stop("Sanity check failure: 0 width variant")
    }
    if (nrow(vcf) > 0 && all(rowData(vcf)@ranges@width == 1)) {
      browser()
      stop("Sanity check failure: all variants are single base insertions")
    }
    if (nrow(vcf) > 0 && all(rowData(vcf)@ranges@width == 1)) {
      browser()
      stop("Sanity check failure: 0 length SV")
    }
    return (vcf)
  }
  vcfs <- lapply(vcfFileList, LoadVcf)
  names(vcfs) <- GetMetadataId(vcfFileList)
  vcfs[sapply(vcfs, is.null)] <- NULL # Remove NULL VCFs list
  write(paste("Loaded", length(vcfs), "VCFs"), stderr())
  return(vcfs)
}

VcfGet <- list()
VcfGet$SVLEN <- function(vcf) {
  if (nrow(vcf) == 0) return(integer(0))
	svcol <- info(vcf)$SVLEN
	if (!is.null(svcol)) {
    return(sapply(svcol, "[", 1))  # grab first allele
	}
  if (attr(vcf, "metadata")$CX_CALLER %in% c("samtools", "gatk", "dindel", "SOAPindel", "soapindel")) {
		# calculate from genotype
    # grab the length of the first alt allele
    # Notes:
    # unlist flattens multiple alleles
    # cumsum of elementLenths returns offset of the last allele in the flattened list
    # to get the first, we shift over by 1 so [1, 3, 4] becomes [+1+, 2, 4, -5-]
    return(width(unlist(alt(vcf))[c(1, 1 + head(cumsum(elementLengths(alt(vcf))), -1))])  - elementLengths(ref(vcf)))
		#return(unlist(lapply(alt(vcf), function (stringSet) { elementLengths(stringSet)[1] } )) - elementLengths(ref(vcf)))
	}
  stop("Unable to determine SVLEN from INFO columns")
}
VcfGet$SVTYPE <- function(vcf) {
  svcol <- info(vcf)$SVTYPE
  if (!is.null(svcol)) {
		return(svcol)
	} else if (length(fixed(vcf)) == 0) {
		return(c())
	} else {
		# fall-back to infer from SVLEN
		return(ifelse(VcfGet$SVLEN(vcf) > 0, "INS", "DEL"))
	}
}
# VCF callers
GenerateDetailedCallInformation <- function(metadata, vcfs) {
	write("Constructing indel caller accuracy data frames", stderr())
	vcfFileList <- list.files(".", pattern="[0-9a-f]{32}.vcf$")
	ProcessVcf <- function (filename) {
    if (!any(metadata$Id == GetMetadataId(filename))) {
      write(paste("Missing metadata for ", filename), stderr())
      return(NULL)
    }
		write(paste("Processing", filename), stderr())
		result <- list()
		result$md <- metadata[GetMetadataId(filename),]
		result$vcf <- vcfs[[result$md$Id]]
		result$referenceVcf <- vcfs[[GetMetadataId(result$md$CX_REFERENCE_VCF)]]
    if (is.null(result$vcf)) {
      return(NULL)
    }
    # filter to vcfs to only those that pass filters
		if (result$md$CX_CALLER %in% c("varscan2")) {
        # TODO: update simulator to improve strandedness of simulated reads
		} else {
		  result$vcf <- result$vcf[rowData(result$vcf)$FILTER %in% c(".", "PASS", "pass"), ]
		}
    
		errorLocationBases <- max(
      result$md$CX_READ_FRAGMENT_LENGTH - 2 * result$md$CX_READ_LENGTH,
      result$md$CX_READ_FRAGMENT_LENGTH / 10) # 1 stddev
		if (result$md$CX_CALLER %in% c(
		  "gatk", "samtools", "soapindel", "SOAPindel", "varscan2", "dindel",
		  "pindel", "splitread", "svseq", "hydra", "socrates")) {
		  # exact breakpoint callers should call exactly
		  # give some margin since callers such as samtools can call strangly.
		  # Eg: CTTA CTT for a 1 base deletion
		  errorLocationBases <- 10
		}
    errorLengthBase <- errorLocationBases
		errorLengthMargin <- 0.20
		refSvType <- VcfGet$SVTYPE(result$referenceVcf)
		refSvLen <- VcfGet$SVLEN(result$referenceVcf)
		callerSvType <- VcfGet$SVTYPE(result$vcf)
		callerSvLen <- VcfGet$SVLEN(result$vcf)
		
		# findOverlaps does not required VCF to be sorted
		
		hits <- findOverlaps(rowData(result$vcf), rowData(result$referenceVcf), maxgap=errorLocationBases, type="equal", ignore.strand=TRUE)
		result$tp <- as.data.frame(hits)
		result$tp <- rename(result$tp, c("queryHits"="vcfIndex", "subjectHits"="refIndex"))
		result$tp$svType <- refSvType[result$tp$refIndex]
		result$tp$svLen <- refSvLen[result$tp$refIndex]
		result$tp$svTypeCalled <- callerSvType[result$tp$vcfIndex]
		result$tp$svLenCalled <- callerSvLen[result$tp$vcfIndex]
		result$tp$svLenDiff <- abs(result$tp$svLen - result$tp$svLenCalled)
		result$tp$qual <- qual(result$vcf)[result$tp$vcfIndex]
		# filter hits that call the wrong SV
		result$tp <- result$tp[
		  sign(result$tp$svLenCalled)==sign(result$tp$svLen) &
		    abs(result$tp$svLenCalled - result$tp$svLen) < errorLengthBase + abs(result$tp$svLen) * errorLengthMargin & 
		    abs(result$tp$svLenCalled) - abs(result$tp$svLen) < abs(result$tp$svLen) # can't be out by more than the size of the SV
		  , ]
	  # filter hits matching multiple SVs to the one closest in length
	  result$tp <- result$tp[order(result$tp$svLenDiff), ]
	  result$tp <- result$tp[!duplicated(result$tp$vcfIndex), ]
    
		# sort so highest quality matches are first
		result$tp <- result$tp[order(result$tp$qual, decreasing=TRUE), ]
		result$tp$isDuplicate <- duplicated(result$tp$refIndex)
		# is there a better way than using %in% to get the variant indexes not in the hits list?
		referenceFalseNegative <- which(!(seq_len(nrow(result$referenceVcf)) %in% result$tp$refIndex))
		callerFalsePositive <- which(!(seq_len(nrow(result$vcf)) %in% result$tp$vcfIndex))
		
		# False Negatives
		result$fn <- data.frame(
		  refIndex=referenceFalseNegative,
		  svType=refSvType[referenceFalseNegative],
		  svLen=refSvLen[referenceFalseNegative],
		  row.names=NULL
		)
		# False Positives
		result$fp <- data.frame(
		  vcfIndex=callerFalsePositive,
		  svTypeCalled=callerSvType[callerFalsePositive],
		  svLenCalled=callerSvLen[callerFalsePositive],
		  qual=qual(result$vcf)[callerFalsePositive],
		  row.names=NULL
		)
		# for now we treat the length of FP as equal to the length of the closest reference SV
		buckets <- sort(unique(abs(unlist(info(result$referenceVcf)$SVLEN))))
		result$fp$svType=result$fp$svTypeCalled
		result$fp$svLen=sign(result$fp$svLenCalled) * buckets[findInterval(abs(result$fp$svLenCalled), (buckets + c(1, head(buckets, -1))) / 2)]
		
		SetCommonFields <- function(df, md) {
		  len <- nrow(df)
		  df$caller <- rep(md$CX_CALLER, len)
		  df$aligner <- rep(md$CX_ALIGNER, len)
		  df$softclipped <- rep(md$CX_ALIGNER_SOFTCLIP == 1, len)
		  df$fragLength <- rep(md$CX_READ_FRAGMENT_LENGTH, len)
		  df$readLength <- rep(md$CX_READ_LENGTH, len)
		  df$readDepth <- rep(md$CX_READ_DEPTH, len)
		  df$id <- rep(md$Id, len)
		  return (df)
		}
		result$tp <- SetCommonFields(result$tp, result$md)
		result$fp <- SetCommonFields(result$fp, result$md)
		result$fn <- SetCommonFields(result$fn, result$md)  	
		# Sanity checks
		if (nrow(result$tp) != length(unique(result$tp$vcfIndex))) {
		  browser()
		  stop("Sanity check failure: call should not match multiple reference SVs. Reference SVs are non-overlapping")
		}
		if (sum(!result$tp$isDuplicate) + nrow(result$fn) != nrow(result$referenceVcf)) {
		  browser()
		  stop("Sanity check failure: TP + FN != reference count")
		}
		if (nrow(result$fp) + nrow(result$tp) != nrow(result$vcf)) {
		  browser()
		  stop("Sanity check failure: FP + TP != VCF SV record count")
		}
		return(result)
	}
	vcfdfs <- lapply(vcfFileList, ProcessVcf)
	vcfdfs[sapply(vcfdfs, is.null)] <- NULL # Remove NULL VCFs list
	resultDfs <- list()
	# aggregate data frames
	resultDfs$fn <- do.call("rbind", lapply(vcfdfs, "[[", "fn"))
	resultDfs$fp <- do.call("rbind", lapply(vcfdfs, "[[", "fp"))
	resultDfs$tp <- do.call("rbind", lapply(vcfdfs, "[[", "tp"))
	return(resultDfs)
}

ApplyDetailAggregation <- function (f, detail, ...) {
	# assume final column is the aggregate column
	result <- f(detail$tp[detail$tp$isDuplicate == FALSE, ], ...)
	stopifnot(!colnames(result)[length(result)] %in% colnames(detail))
	colnames(result)[length(result)] <- "tp"
	result <- merge(result, f(detail$tp[detail$tp$isDuplicate == TRUE, ], ...), all=TRUE)
	stopifnot(!colnames(result)[length(result)] %in% colnames(detail))
	colnames(result)[length(result)] <- "tpDup"
	result <- merge(result, f(detail$fn, ...), all=TRUE)
	stopifnot(!colnames(result)[length(result)] %in% colnames(detail))
	colnames(result)[length(result)] <- "fn"
	result <- merge(result, f(detail$fp, ...), all=TRUE)
	stopifnot(!colnames(result)[length(result)] %in% colnames(detail))
	colnames(result)[length(result)] <- "fp"
	# NA = 0
	result$tp[is.na(result$tp)] <- 0
	result$tpDup[is.na(result$tpDup)] <- 0
	result$fn[is.na(result$fn)] <- 0
	result$fp[is.na(result$fp)] <- 0
	# we have an unreasonable number of true negatives
	# TODO(cameron.d): are we handling TP duplicates correctly?
	result$reftotal <- result$tp + result$fn
	result$sensitivity <- result$tp / (result$tp + result$fn)
	#result$fpr <- result$fp / (result$fp + result$tn)
	#result$accuracy <- (result$tp + result$tn) / (result$tp + result$tn + result$fp + result$fn)
	#result$specificity <- result$fp / (result$fp + result$tn)
	result$precision <- result$tp / (result$tp + result$tpDup + result$fp)
	#result$negativePredictiveValue <- result$tn / (result$tn + result$fn)
	result$fdr <- (result$fp + result$tpDup) / (result$fp + result$tp)
	return(result)
}

LoadTruth <- function(vcf, truthvcf, annotationFunction = function(vcf) info(vcf), ...) {
  if (any(rowData(vcf)$QUAL <= 0)) {
    stop("Precondition failure: not all variants have positive quality score")
  }
  vcfgr <- vcftobpgr(vcf)
  vcfdf <- annotationFunction(vcf)
  truthgr <- vcftobpgr(truthvcf)
  truthdf <- as.data.frame(truthgr)
  hits <- breakpointHits(vcfgr, truthgr, ...)
  hits$qual <- rowData(vcf)$QUAL[hits$queryHits]
  hits <- hits[order(hits$qual),]
  # set up true
  vcfdf$tp <- FALSE
  vcfdf$tp[hits$queryHits] <- TRUE
  truthdf$tp <- FALSE
  truthdf$tp[hits$subjectHits] <- TRUE
  truthdf$qual <- 0
  truthdf$qual[hits$subjectHits] <- hits$qual
  #vcfdf <- cbind(vcfdf, vcf@metadata)
  #truthdf <- cbind(truthdf, vcf@metadata)
  return(list(vcf=vcf, calls=vcfdf, truth=truthdf))
}
LoadTruthSummary <- function(vcfs, ...) {
  truthSet <- lapply(vcfs, function(vcf) {
    md <- attr(vcf, "metadata")
    if (is.null(md)) return(NULL)
    if (is.null(md$CX_CALLER)) return(NULL)
    truthvcf <- vcfs[[GetMetadataId(md$CX_REFERENCE_VCF)]]
    if (is.null(truthvcf)) {
      warning(paste0("Missing truth for", md$Id))
    }
    write(paste0("Loading truth for ", md$Id), stderr())
    callTruthPair <- LoadTruth(vcf, truthvcf, ...)
    callTruthPair$calls <- cbind(callTruthPair$calls, md)
    callTruthPair$truth <- cbind(callTruthPair$truth, md)
    return(callTruthPair)
  })
  truthSet[sapply(truthSet, is.null)] <- NULL # Remove NULLs
  calls = rbindlist(sapply(truthSet, function(x) x$calls), use.names=TRUE, fill=TRUE)
  truth = rbindlist(sapply(truthSet, function(x) x$truth), use.names=TRUE, fill=TRUE)
  return(list(calls=calls, truth=truth))
}

