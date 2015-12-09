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
theme_set(theme_bw())

scale_y_power4 <- scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1)^4, labels=c("0", "", "", "", "", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"))
scale_y_power5 <- scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1)^5, labels=c("0", "", "", "", "", "0.5", "", "0.7", "0.8", "0.9", "1.0"))

# filters the data frame to the subset of calls to be used in the body of the gridss paper
filter_main_callers <- function(df, col="CX_CALLER") {
  df <- df[!(str_detect(df[[col]], "gridss[/].*[/]")),]
  df[[col]] <- factor(str_extract(df[[col]], "^[^/]+"))
  df <- df[df[[col]] %in% c("gridss", "breakdancer", "delly", "lumpy", "pindel", "socrates"),]
  df[[col]] <- relevel(df[[col]], "gridss")
  return(df)
}

saveplot <- function(file=file, ...) {
  ggsave(paste0("png/", file, ".png"), ...)
  ggsave(paste0("eps/", file, ".eps"), ...)
}


findOverlaps_type_equal_df <- function(query, subject, ...) {
  shits <- findOverlaps(query, subject, type="start", ...)
  ehits <- findOverlaps(query, subject, type="end", ...)
  hits <- data.frame(queryHits=c(queryHits(shits), queryHits(ehits)), subjectHits=c(subjectHits(shits), subjectHits(ehits)))
  hits <- hits[duplicated(hits),]
}