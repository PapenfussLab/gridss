library(VariantAnnotation)
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
gridss.truthdetails.processvcf.vcftodf <- function(vcf) {
  i <- info(vcf)
  df <- data.frame(variantid=names(rowData(vcf)))
  df$QUAL <- fixed(vcf)$QUAL
  df$EVENT <-i$EVENT
  df$SOMATIC <- i$SOMATIC
  matchLength <- as.integer(gsub(".*_", "", as.character(i$TRUTH_MATCHES)))
  mismatchLength <- as.integer(gsub(".*_", "", as.character(i$TRUTH_MISREALIGN)))
  df$SVLEN <- ifelse(is.na(matchLength), mismatchLength, matchLength)
  df$SVTYPE <- i$SVTYPE
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
}



