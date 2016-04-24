#
# Checks to see which callers have/have not been run for which inputs
#
dataset <- "fastcompare"

setwd(paste0(ifelse(as.character(Sys.info())[1] == "Windows", "W:/", "~/"), "dev/gridss/src/test/r"))
source("common.R")

setwd(paste0(ifelse(as.character(Sys.info())[1] == "Windows", "W:/", "~/"), "i/data.", dataset))
metadata <- data.table(LoadMetadata(), key="Id")
metadata <- metadata[!is.na(metadata$CX_CALLER),]
vcfdf <- data.table(vcf=list.files(".", "*.vcf$"))
vcfdf$Id <- GetMetadataId(vcfdf$vcf)
vcfdf$vcfsize <- file.size(vcfdf$vcf)
setkey(vcfdf, "Id")
metadata <- vcfdf[metadata]
metadata$CX_ALIGNER <- ifelse(is.na(metadata$CX_ALIGNER), ".", as.character(metadata$CX_ALIGNER))
metadata$lock <- sapply(metadata$Id, function (id) {
  filename <- paste0(id, ".lock")
  if (file.exists(filename)) return(readChar(filename, file.info(filename)$size - 1))
  return("")
})

sliceby <- c("CX_ALIGNER", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_FRAGMENT_STDDEV",  "CX_READ_LENGTH", "CX_REFERENCE_VCF_VARIANTS")
slice <- unique(metadata[, sliceby, with=FALSE])
output <- metadata[, merge(slice, .SD, by=sliceby, all=TRUE), by="CX_CALLER"]
output$caller <- str_extract(output$CX_CALLER, "^[^/]+")
# cortex runs off FASTA files (via bwamem output so downsampling uses the same reads as the bwamem downsampling)
output <- output[!(output$caller=="cortex" & output$CX_ALIGNER != "bwamem"),]
output <- output[!(output$caller !="variationhunter" & output$CX_ALIGNER == "mrfast/2.6.0.1"),]
output <- output[!(output$caller =="variationhunter" & output$CX_ALIGNER != "mrfast/2.6.0.1"),]
output <- output[!(!(output$caller %in% c("crest", "socrates")) & output$CX_ALIGNER == "bowtie2"),]
write.csv(output, paste0(ifelse(as.character(Sys.info())[1] == "Windows", "W:/", "~/"), "i/progress.", dataset, ".csv"))

