#!/usr/bin/env Rscript
library(argparser)
argp = arg_parser("Filters a raw GRIDSS VCF into somatic call subsets.")
argp = add_argument(argp, "--pondir", default=NA, help="Directory containing Panel Of Normal bed/bedpe used to filter FP somatic events. USer create_gridss_pon.R to generate the PON.")
argp = add_argument(argp, "--ref", default="BSgenome.Hsapiens.UCSC.hg19", help="Reference genome to use. Must be a valid installed BSgenome package")
argp = add_argument(argp, "--input", help="GRIDSS VCF")
argp = add_argument(argp, "--output", help="High confidence somatic subset")
argp = add_argument(argp, "--fulloutput", help="Full call set excluding obviously germline call.")
argp = add_argument(argp, "--normalordinal", type="integer", default=1, help="Ordinal of matching normal sample in the VCF")
argp = add_argument(argp, "--scriptdir", default=ifelse(sys.nframe() == 0, "./", dirname(sys.frame(1)$ofile)), help="Path to libgridss.R script")
argp = add_argument(argp, "--gc", flag=TRUE, help="Perform garbage collection after freeing of large objects. ")
# argv = parse_args(argp, argv=c("--input", "D:/hartwig/down/COLO829R_COLO829T.gridss.vcf", "--output", "D:/hartwig/temp/out.vcf", "-f", "D:/hartwig/temp/full.vcf", "-p", "D:/hartwig/dbs/gridss/pon3792v1", "--scriptdir", "D:/hartwig/scripts/gridss", "--gc"))
#argv = parse_args(argp, argv=c("--input", "D:/hartwig/down/annotate_variants.vcf", "--output", "D:/hartwig/temp/out.vcf", "-f", "D:/hartwig/temp/full.vcf", "-p", "D:/hartwig/dbs/gridss/pon3792v1", "--scriptdir", "D:/hartwig/scripts/gridss", "--gc"))
argv = parse_args(argp)

if (!file.exists(argv$input)) {
  msg = paste(argv$input, "not found")
  write(msg, stderr())
  print(argp)
  stop(msg)
}
if (is.na(argv$pondir)) {
  argv$pondir = NULL
} else if (!dir.exists(argv$pondir)) {
  msg = paste(argv$pondir, "not found")
  write(msg, stderr())
  print(argp)
  stop(msg)
}
refgenome=eval(parse(text=paste0("library(", argv$ref, ")\n", argv$ref)))

library(tidyverse)
library(readr)
library(stringr)
libgridssfile = paste0(argv$scriptdir, "/", "libgridss.R")
if (file.exists(libgridssfile)) {
  tmpwd = getwd()
  setwd(argv$scriptdir)
  source("libgridss.R")
  setwd(tmpwd)
} else {
  msg = paste("Could not find libgridss.R in", argv$scriptdir, " - please specify a --scriptdir path to a directory containing the required scripts")
  write(msg, stderr())
  print(argp)
  stop(msg)
}


# Filter to somatic calls
write(paste(Sys.time(), "Reading", argv$input), stderr())
raw_vcf = readVcf(argv$input)
tumourordinal = seq(ncol(geno(raw_vcf)$VF))[-argv$normalordinal]
# hard filter variants that are obviously not somatic
full_vcf = raw_vcf[geno(raw_vcf)$QUAL[,argv$normalordinal] / VariantAnnotation::fixed(raw_vcf)$QUAL < 4 * gridss.allowable_normal_contamination]
rm(raw_vcf)
if (argv$gc) { gc() }
# hard filter unpaired breakpoints (caused by inconsistent scoring across the two breakends)
full_vcf = full_vcf[is.na(info(full_vcf)$PARID) | info(full_vcf)$PARID %in% names(full_vcf)]
full_vcf = align_breakpoints(full_vcf)
# Add header fields
full_vcf = addVCFHeaders(full_vcf)

info(full_vcf)$BPI_AF = rep("", length(full_vcf))
filters = rep("", length(full_vcf))
names(filters) = names(full_vcf)

write(paste(Sys.time(), "Parsing single breakends", argv$input), stderr())
begr = breakpointRanges(full_vcf, unpartneredBreakends=TRUE)
write(paste(Sys.time(), "Calculating single breakend VAF", argv$input), stderr())
begr$af = round(gridss_be_af(begr, full_vcf, tumourordinal), 5)
begr$af_str = as.character(begr$af)
if (length(begr) > 0) {
  info(full_vcf[names(begr)])$BPI_AF = begr$af_str
}
write(paste(Sys.time(), "Filtering single breakends", argv$input), stderr())
befiltered = gridss_breakend_filter(begr, full_vcf, pon_dir=argv$pondir, normalOrdinal=argv$normalordinal, tumourOrdinal=tumourordinal)
filters[names(begr)] = befiltered
rm(befiltered)
full_vcf = full_vcf[passes_very_hard_filters(filters)]
filters = filters[passes_very_hard_filters(filters)]
begr = begr[names(begr) %in% names(full_vcf)]
if (argv$gc) { gc() }


write(paste(Sys.time(), "Parsing breakpoints", argv$input), stderr())
bpgr = breakpointRanges(full_vcf, unpartneredBreakends=FALSE)
write(paste(Sys.time(), "Calculating breakpoint VAF", argv$input), stderr())
bpgr$af = round(gridss_bp_af(bpgr, full_vcf, tumourordinal), 5)
bpgr$af_str = paste(bpgr$af, partner(bpgr)$af, sep=",")
if (length(bpgr) > 0) {
  info(full_vcf[names(bpgr)])$BPI_AF = bpgr$af_str
}
write(paste(Sys.time(), "Filtering breakpoints", argv$input), stderr())
bpfiltered = gridss_breakpoint_filter(bpgr, full_vcf, bsgenome=refgenome, pon_dir=argv$pondir, normalOrdinal=argv$normalordinal, tumourOrdinal=tumourordinal)
filters[names(bpgr)] = bpfiltered
if (argv$gc) { gc() }

# shadow breakpoint removed due to initial mapq20 filter reducing FP rate
# bpfiltered = .addFilter(bpfiltered, "shadow", is_shadow_breakpoint(bpgr, begr, full_vcf))

#bpfiltered = .addFilter(bpfiltered, "LOW_Qual", bpgr$QUAL < gridss.min_qual)
#som_llr = gridss_breakpoint_somatic_llr(full_vcf, bpgr=bpgr, contamination_rate=gridss.allowable_normal_contamination)

# - filter to only decent length assemblies?
#begr$calls_1k_window = countOverlaps(begr, rowRanges(full_vcf), ignore.strand=TRUE, maxgap=1000)

# Remove very hard filtered variants
full_vcf = full_vcf[passes_very_hard_filters(filters)]
unpaired_breakpoint = !is.na(info(full_vcf)$PARID) & !(info(full_vcf)$PARID %in% names(full_vcf))
full_vcf = full_vcf[!unpaired_breakpoint]
filters = filters[passes_very_hard_filters(filters)]
filters = filters[!unpaired_breakpoint]
rm(unpaired_breakpoint)

vcf = full_vcf[passes_soft_filters(filters)]
vcf = vcf[is.na(info(vcf)$PARID) | info(vcf)$PARID %in% names(vcf)]
bpgr = bpgr[names(bpgr) %in% names(vcf)]
if (argv$gc) { gc() }

write(paste(Sys.time(),"Calculating transitive links", argv$input), stderr())
# transitive calling
transitive_df = transitive_calls(vcf, bpgr, report="max2") %>%
  # only make transitive calls were we actually know the path
  filter(!has_multiple_paths) %>%
  mutate(type="transitive")
# now we filter imprecise variants
is_imprecise = !(is.na(info(vcf)$IMPRECISE) | !info(vcf)$IMPRECISE) |
  !((!is.na(info(vcf)$PARID) & info(vcf)$ASSR + info(vcf)$SR + info(vcf)$IC > 0) |
      (is.na(info(vcf)$PARID) & info(vcf)$BASSR + info(vcf)$BSC > 0))
filters[names(vcf)[is_imprecise]] = paste0(filters[names(vcf)[is_imprecise]], ";imprecise")
filters[transitive_df$linked_by] = paste0(filters[transitive_df$linked_by], ";transitive")

vcf = full_vcf[passes_soft_filters(filters)]
vcf = vcf[is.na(info(vcf)$PARID) | info(vcf)$PARID %in% names(vcf)]
bpgr = bpgr[names(bpgr) %in% names(vcf)]
begr = begr[names(begr) %in% names(vcf)]
if (argv$gc) { gc() }

asm_linked_df = NULL
write(paste(Sys.time(),"Calculating assembly links", argv$input), stderr())
# Assembly-based event linking
if (length(vcf) > 0) {
  asm_linked_df = linked_assemblies(vcf) %>%
    mutate(type="asm")
}

link_df = bind_rows(asm_linked_df, transitive_df) %>%
  mutate(linking_group=str_replace(linked_by, "/.*$", "")) %>%
  mutate(pass=passes_final_filters(vcf[vcfId])) %>%
  group_by(linking_group) %>%
  mutate(pass=any(pass)) %>%
  ungroup() %>%
  filter(pass)

write(paste(Sys.time(),"Calculating bebe insertion links", argv$input), stderr())
# Insertion linkage
bebeins_link_df = linked_by_breakend_breakend_insertion_classification(begr) %>%
  group_by(linked_by) %>%
  mutate(pass=passes_final_filters(vcf[vcfId])) %>%
  mutate(pass=any(pass)) %>%
  ungroup() %>%
  filter(pass) %>%
  mutate(type="bebeins")
write(paste(Sys.time(),"Calculating bebp insertion links", argv$input), stderr())
bebpins_link_df = linked_by_breakpoint_breakend_insertion_classification(bpgr, begr) %>%
  group_by(linked_by) %>%
  mutate(pass=passes_final_filters(vcf[vcfId])) %>%
  mutate(pass=any(pass)) %>%
  ungroup() %>%
  filter(pass) %>%
  mutate(type="bebpins")
# Inversion linkage
write(paste(Sys.time(),"Calculating simple inversions", argv$input), stderr())
inv_link_df = linked_by_simple_inversion_classification(bpgr) %>%
  group_by(linked_by) %>%
  mutate(pass=passes_final_filters(vcf[vcfId])) %>%
  mutate(pass=any(pass)) %>%
  ungroup() %>%
  filter(pass) %>%
  mutate(type="inv")
# Deletion bridge linkage
# TODO: do we want to do this?
# I'm suspicious of the model used in ChainFinder PMC3690918
# Notably: I'm suspicous that "repair with major DNA loss" is actually a thing
# given the catastrophic nature of chromo*, a more reasonable explaination is
# an additional DSB with the subsequent loss of that DNA fragment.
# Given the focal nature of chromoplexy, ChainFinder works because it just
# finds the focal events, not because the model is correct.
# TODO: show this by modelling additional focal DSBs
write(paste(Sys.time(),"Calculating dsb links", argv$input), stderr())
dsb_link_df = linked_by_dsb(bpgr) %>%
  group_by(linked_by) %>%
  mutate(pass=passes_final_filters(vcf[vcfId])) %>%
  mutate(pass=any(pass)) %>%
  ungroup() %>%
  filter(pass) %>%
  mutate(type="dsb")

write(paste(Sys.time(),"Removing duplicated/conflicting links", argv$input), stderr())
# linking priorities:
# - asm independent of other linkages
# - transitive independent of other linkages
# - ins, inv, dsb linkages
event_link_df = bind_rows(
  bebeins_link_df,
  bebpins_link_df,
  inv_link_df,
  dsb_link_df) %>%
  dplyr::select(vcfId, linked_by) %>%
  mutate(
    QUAL=rowRanges(vcf)[vcfId]$QUAL,
    hasPolyA=str_detect(rowRanges(vcf[vcfId])$ALT, "A{16}")) %>%
  group_by(linked_by) %>%
  # filter events where supporting fragment counts differ by too much
  mutate(
    max_supporting_fragment_count = max(ifelse(is.na(info(full_vcf[vcfId])$PARID), info(full_vcf[vcfId])$BVF, info(full_vcf[vcfId])$VF)),
    min_supporting_fragment_count = min(ifelse(is.na(info(full_vcf[vcfId])$PARID), info(full_vcf[vcfId])$BVF, info(full_vcf[vcfId])$VF)),
    hasPolyA=any(hasPolyA)
    ) %>%
  filter(min_supporting_fragment_count >= gridss.min_rescue_portion * max_supporting_fragment_count | hasPolyA)

write(paste(Sys.time(),"Calculating final linkage annotation", argv$input), stderr())
# Only keep the best QUAL event linkage
event_link_df = event_link_df %>%
  group_by(linked_by) %>%
  mutate(linkQUAL = pmin(QUAL)) %>%
  group_by(vcfId) %>%
  filter(QUAL == linkQUAL) %>%
  group_by(linked_by)
# Don't event link to PON filtered variants
event_link_df = event_link_df %>%
  filter(!str_detect(filters[vcfId], "PON"))
# Fix up pairing
event_link_df = event_link_df %>%
  filter(n() == 2) %>%
  ungroup()

# include both breakends of any linked breakpoints
# as linkage can be breakend specific (e.g. assembly, bpbeins)
link_rescue = bind_rows(link_df, event_link_df) %>% pull(vcfId) %>% unique()
link_rescue = c(link_rescue, bpgr[link_rescue[link_rescue %in% names(bpgr)]]$partner)

# Note that we don't rescue equivalent events
begr$partner = rep(NA, length(begr))
eqv_link_df = linked_by_equivalent_variants(full_vcf, as(rbind(as.data.frame(bpgr), as.data.frame(begr)), "GRanges"), bsgenome=refgenome) %>%
  filter(passes_final_filters(vcf[vcfId]) | vcfId %in% link_rescue) %>%
  group_by(linked_by) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  mutate(type="eqv")

link_summary_df = bind_rows(link_df, event_link_df, eqv_link_df) %>%
  group_by(vcfId) %>%
  summarise(linked_by=paste0(linked_by, collapse=","))

# Add linking information
info(full_vcf)$LOCAL_LINKED_BY = rep("", length(full_vcf))
info(full_vcf)$REMOTE_LINKED_BY = rep("", length(full_vcf))
info(full_vcf[link_summary_df$vcfId])$LOCAL_LINKED_BY = link_summary_df$linked_by
info(full_vcf[!is.na(info(full_vcf)$PARID)])$REMOTE_LINKED_BY = info(full_vcf[info(full_vcf[!is.na(info(full_vcf)$PARID)])$PARID])$LOCAL_LINKED_BY

write(paste(Sys.time(),"Final qual filtering ", argv$output), stderr())
# final qual filtering
fails_qual_without_rescue = !passes_final_filters(full_vcf, include.existing.filters=FALSE) & !(names(full_vcf) %in% link_rescue)
filters[names(full_vcf)[fails_qual_without_rescue]] = paste0(filters[names(full_vcf)[fails_qual_without_rescue]], ";qual")

################
# Write outputs
VariantAnnotation::fixed(full_vcf)$FILTER = ifelse(str_remove(filters, "^;") == "", "PASS", str_remove(filters, "^;"))
if (argv$gc) { gc() }
if (!is.na(argv$output)) {
  write(paste(Sys.time(),"Writing ", argv$output), stderr())
  vcf = full_vcf[passes_soft_filters(filters)]
  vcf = vcf[is.na(info(vcf)$PARID) | info(vcf)$PARID %in% names(vcf)]
  writeVcf(vcf, argv$output, index=TRUE)
}
if (!is.na(argv$fulloutput)) {
  write(paste(Sys.time(),"Writing ", argv$fulloutput), stderr())
  vcf = full_vcf[passes_very_hard_filters(filters)]
  vcf = vcf[is.na(info(vcf)$PARID) | info(vcf)$PARID %in% names(vcf)]
  writeVcf(vcf, argv$fulloutput, index=TRUE)
}





