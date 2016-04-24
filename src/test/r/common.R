##################
# Processing steps (64 bit linux required)
##################
# - Create environment modules that put all executables in PATH for the following tools:
# -- bwa/0.7.12
# -- lumpy/0.2.11
# -- samblaster/0.1.22 (required by lumpy)
# -- ucsc-tools/322
# -- bam-readcount/0.7.4 (required by varscan) [optional]
# -- cap3/20151002 (required by crest) [optional]
# -- crest/0.0.1 [optional]
# -- pindel/0.2.5b6
# -- varscan/2.4.0 [optional]
# -- delly/0.6.8
# -- samtools/1.2
# -- breakdancer/1.4.5
# -- clever/2.0rc3 [optional]
# -- gasv/20140228 [optional]
# Ensure the following tools are on PATH by default:
# -- bwa
# -- bowtie2
# -- novoalign
# -- novosort
# -- samtools
# -- ucsc-tools
# -- Picard tools (with wrappers for DownsampleSam and SortSam)
# - Download Socartes 1.12 to ~/src/socrates/1.12/target/socrates-1.12-jar-with-dependencies.jar
# - Download gridss to ~/bin/gridss-0.9.2-jar-with-dependencies.jar
# - Download hg19.fa to ~/reference_genomes/human/hg19.fa
# --  create hg19.fa.fai, hg19.dict, and bwa, bowtie2, novoalign indexes using hg19.fa as the prefix, and a BLAT index
# - Download UCSC repeatmaster track and decompress so that ~/Papenfuss_lab/projects/reference_genomes/human/hg19/UCSC/chromOut/12/chr12.fa.out exists
# - Create symbolic links 
# - git clone http://github.com/Papenfuss_Lab/gridss
# - cd src/test/sim
# - create symbolic link from src/test/sim to ~/i
# - Download art 1.51 to ~/i/tools/art
library(ggplot2)
library(rtracklayer)

theme_set(theme_bw())

rootdir <- ifelse(as.character(Sys.info())[1] == "Windows", "W:/", "~/")

setwd(paste0(rootdir, "dev/gridss/src/test/r"))
source("libvcf.R")

scale_y_power4 <- scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1)^4, labels=c("0", "", "", "", "", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"))
scale_y_power5 <- scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1)^5, labels=c("0", "", "", "", "", "0.5", "", "0.7", "0.8", "0.9", "1.0"))

# filters the data frame to the subset of calls to be used in the body of the gridss paper
main_callers_subset <- function(metadata,
                                slice=c("CX_CALLER", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_FRAGMENT_STDDEV",	"CX_READ_LENGTH", "CX_REFERENCE_VCF"),
                                callers=c("gridss", "breakdancer", "delly", "lumpy", "pindel", "socrates")) { #, "crest", "cortex"
  slice <- slice[slice %in% names(metadata)] # only pull out columns that actually exist
  #df <- df[!(str_detect(df[[col]], "gridss[/].*[/]")),]
  metadata$CX_CALLER <- str_extract(metadata$CX_CALLER, "^[^/]+")
  metadata <- metadata[is.na(metadata$CX_CALLER) | metadata$CX_CALLER %in% callers,]
  if ("CX_BAM" %in% names(metadata)) {
    # Socrates uses the bowtie2 version if it exists regardless of initial aligner
    metadata$CX_CALLER[!is.na(metadata$CX_CALLER) & !is.na(metadata$CX_BAM) & str_detect(metadata$CX_BAM, ".bam.bt2.bam") & metadata$CX_CALLER == "socrates"] <- "bowtie2"
  }
  if ("CX_ALIGNER" %in% names(metadata)) {
    callermd <- metadata[!is.na(metadata$CX_CALLER),]
    callermd$preferredcaller <- 
      callermd$CX_CALLER %in% c("cortex") |
      (callermd$CX_ALIGNER == "bwamem" & callermd$CX_CALLER %in% c("gridss", "breakdancer", "delly", "lumpy", "pindel", "gasv", "varscan", "clever", "crest")) | 
      (callermd$CX_ALIGNER == "bowtie2" & callermd$CX_CALLER %in% c("socrates"))    
    mdpreferredstatus <- data.table(callermd)[, list(haspreferred=any(preferredcaller)), by=slice]
    callermd <- merge(callermd, data.frame(mdpreferredstatus), by=slice, suffixes=c("", ""))
    callermd <- callermd[(callermd$preferredcaller & callermd$haspreferred) | !callermd$haspreferred,]
    if (any(!callermd$haspreferred)) {
      warning(paste("Using fallback aligner for ", paste(unique(callermd$CX_CALLER[!callermd$haspreferred]), collapse=" ")))
    }
    callermd$preferredcaller <- NULL
    callermd$haspreferred <- NULL
    metadata <- rbind(metadata[is.na(metadata$CX_CALLER),], callermd)
  }
  metadata$CX_CALLER <- as.factor.gridss.first(metadata$CX_CALLER)
  if (anyDuplicated(metadata[!is.na(metadata$CX_CALLER), slice])) {
    md <- metadata[!is.na(metadata$CX_CALLER), slice]
    md[duplicated(md),]
    browser()
    stop("Found multiple metadata instance for main caller")
  }
  return(metadata)
}

text_default_calls <- "Default calls"
text_all_calls <- "All calls"
prettyFilter <- function(x) factor(x, levels=c(text_default_calls, text_all_calls))

as.factor.gridss.first <- function(x) {
  return(relevel(factor(x), "gridss"))
}

saveplot <- function(file=file, ...) {
  ggsave(paste0("png/", file, ".png"), ...)
  ggsave(paste0("eps/", file, ".eps"), ...)
}

# VariantAnnotation workaround
findOverlaps_type_equal_df <- function(query, subject, ...) {
  shits <- findOverlaps(query, subject, type="start", ...)
  ehits <- findOverlaps(query, subject, type="end", ...)
  hits <- data.frame(queryHits=c(queryHits(shits), queryHits(ehits)), subjectHits=c(subjectHits(shits), subjectHits(ehits)))
  hits <- hits[duplicated(hits),]
}


#
# Metadata and VCF processing
#
GetMetadataId <- function(filenames) {
  return(str_match(as.character(filenames), "[0-9a-f]{32}")[,1])
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
CalculateTruthSummary <- function(vcfs, blacklist=NULL, maxerrorbp, ignoreFilters=TRUE, ...) {
  truthSet <- lapply(vcfs, function(vcf) {
    md <- attr(vcf, "sourceMetadata")
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
    callTruthPair <- CalculateTruth(vcf, truthvcf, blacklist, maxerrorbp, ignoreFilters, ...)
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
# Remove callers without a non-zero output VCF from the metadata listing
WithVcf <- function(metadata, directory=".", pattern="*.vcf$") {
  fileList <- list.files(directory, pattern=pattern)
  zeroSizeFiles = file.info(fileList)$size == 0
  if (any(zeroSizeFiles)) {
    fileList <- fileList[!zeroSizeFiles]
  }
  toRemove <- !(metadata$Id %in% GetMetadataId(fileList))
  if (any(toRemove)) {
    warning(paste("Ignoring", paste(metadata$Id[toRemove], collapse=" "), "due to missing VCFs."))
    return(metadata[!toRemove,])
  }
  return(metadata)
}
# Load VCFs into a list
LoadVcfs <- function(metadata, directory=".", pattern="*.vcf$", existingVcfs=NULL) {
  write("Loading VCFs", stderr())
  fileList <- list.files(directory, pattern=pattern)
  zeroSizeFiles = file.info(fileList)$size == 0
  if (any(zeroSizeFiles)) {
    write(paste("Skipping file", fileList[zeroSizeFiles], "due to 0 size.\n"))
    warning(paste("Skipping files", paste(fileList[zeroSizeFiles]), "due to 0 size.\n"))
    fileList <- fileList[!zeroSizeFiles]
  }
  # exclude already loaded VCFs
  fileList <- fileList[!(GetMetadataId(fileList) %in% names(existingVcfs))]
  # only load VCFS that have metadata
  fileList <- fileList[!is.na(GetMetadataId(fileList)) & GetMetadataId(fileList) %in% metadata$Id]
  #vcfs <- foreach (filename=fileList, .packages="VariantAnnotation") %dopar% { # Parallel load of VCFs
  vcfs <- lapply(fileList, function(filename) {
    id <- GetMetadataId(filename)
    write(paste0("Loading ", filename), stderr())
    vcf <- readVcf(filename, "unknown")
    attr(vcf, "filename") <- filename
    attr(vcf, "id") <- id
    attr(vcf, "sourceMetadata") <- metadata[metadata$Id==id,]
    if (nrow(metadata[metadata$Id==id,]) != 1) {
      browser()
      stop(paste("Metadata error for ", filename))
    }
    return(vcf)
  })
  vcfs[sapply(vcfs, is.null)] <- NULL
  names(vcfs) <- sapply(vcfs, function(vcf) attr(vcf, "id"))
  write(paste("Loaded", length(vcfs), "VCFs"), stderr())
  return(c(existingVcfs, vcfs))
}
