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
  df$somatic <- i$SOMATIC
  df$somaticp <- i$SPV
  matchLength <- as.integer(gsub(".*_", "", as.character(i$TRUTH_MATCHES)))
  mismatchLength <- as.integer(gsub(".*_", "", as.character(i$TRUTH_MISREALIGN)))
  df$svlen <- ifelse(is.na(matchLength), mismatchLength, matchLength)
  df$svtype <- i$SVTYPE
  df$CQUAL <- i$CQUAL
  df$call <- ifelse(!is.na(matchLength), "good", ifelse(!is.na(mismatchLength), "misaligned", "bad"))
  df$LR <- i$LR
  df$LRBP <- i$LRBP
  df <- gridss.truthdetails.processvcf.addtn(df, "RC", i$RC)
  df <- gridss.truthdetails.processvcf.addtn(df, "PC", i$PC)
  df <- gridss.truthdetails.processvcf.addtn(df, "A_BCT", i$A_BCT)
  df$A_BLLM <- i$A_BLLM
  df$A_BLRM <- i$A_BLRM
  df$assemblyLength <- df$A_BLLM + df$A_BLRM
  df$A_EC <- i$A_EC
  df$A_LR <- i$A_LR
  df$A_MQM <- i$A_MQM
  df$A_MQT <- i$A_MQT
  df$A_RM <-i$A_RM
  df <- gridss.truthdetails.processvcf.addtn(df, "A_RP", i$A_RP)
  df <- gridss.truthdetails.processvcf.addtn(df, "A_RPBLM", i$A_RPBLM, max=TRUE)
  df <- gridss.truthdetails.processvcf.addtn(df, "A_SC", i$A_SC)
  df <- gridss.truthdetails.processvcf.addtn(df, "A_SCCLM", i$A_SCCLM)
  df <- gridss.truthdetails.processvcf.addtn(df, "A_SCCLT", i$A_SCCLT)  
  df <- gridss.truthdetails.processvcf.addtn(df, "RPEC", i$RPEC)
  df <- gridss.truthdetails.processvcf.addtn(df, "RPLR", i$RPLR)
  df <- gridss.truthdetails.processvcf.addtn(df, "RPMQLM", i$RPMQLM, max=TRUE)
  df <- gridss.truthdetails.processvcf.addtn(df, "RPMQLT", i$RPMQLT)
  df <- gridss.truthdetails.processvcf.addtn(df, "RPMQRM", i$RPMQRM, max=TRUE)
  df <- gridss.truthdetails.processvcf.addtn(df, "RPMQRT", i$RPMQRT)
  df <- gridss.truthdetails.processvcf.addtn(df, "RPRM", i$RPRM)
  df <- gridss.truthdetails.processvcf.addtn(df, "SCBLRM", i$SCBLRM, max=TRUE)
  df <- gridss.truthdetails.processvcf.addtn(df, "SCBLRT", i$SCBLRT)
  df <- gridss.truthdetails.processvcf.addtn(df, "SCEC", i$SCEC)
  df <- gridss.truthdetails.processvcf.addtn(df, "SCLR", i$SCLR)
  df <- gridss.truthdetails.processvcf.addtn(df, "SCMQRM", i$SCMQRM, max=TRUE)
  df <- gridss.truthdetails.processvcf.addtn(df, "SCMQRT", i$SCMQRT)
  df <- gridss.truthdetails.processvcf.addtn(df, "SCRM", i$SCRM)
}