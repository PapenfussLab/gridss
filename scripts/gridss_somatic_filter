#!/usr/bin/env Rscript
library(argparser)
thisFile <- function() { # https://stackoverflow.com/questions/1815606/determine-path-of-the-executing-script
	cmdArgs <- commandArgs(trailingOnly = FALSE)
	needle <- "--file="
	match <- grep(needle, cmdArgs)
	if (length(match) > 0) {
		# Rscript
		return(normalizePath(sub(needle, "", cmdArgs[match])))
	} else if (is.null(sys.frames()[[1]]$ofile)) {
		return("./")
	} else {
		# 'source'd via R console
		return(normalizePath(sys.frames()[[1]]$ofile))
	}
}
argp = arg_parser("Filters a raw GRIDSS VCF into somatic call subsets.")
argp = add_argument(argp, "--pondir", default=NA, help="Directory containing Panel Of Normal bed/bedpe used to filter FP somatic events. Use gridss.GeneratePonBedpe to generate the PON.")
argp = add_argument(argp, "--ref", default="", help="Reference genome to use. Must be a valid installed BSgenome package")
argp = add_argument(argp, "--input", help="GRIDSS VCF")
argp = add_argument(argp, "--output", help="High confidence somatic subset")
argp = add_argument(argp, "--fulloutput", help="Full call set excluding obviously germline call.")
argp = add_argument(argp, "--plotdir", default="", help="Output directory for plots")
argp = add_argument(argp, "--normalordinal", type="integer", default=1, help="Ordinal of matching normal sample in the VCF")
argp = add_argument(argp, "--tumourordinal", type="integer", nargs=Inf, help="Ordinal of tumour sample(s) in the VCF. Defaults to all samples not listed as matched normals")
argp = add_argument(argp, "--scriptdir", default=thisFile(), help="Path to libgridss.R script")
argp = add_argument(argp, "--configdir", default=".", help="Path to gridss.config.R script relative to scriptdir. Defaults to '.' (same directory as libgridss.R)")
argp = add_argument(argp, "--gc", flag=TRUE, help="Perform garbage collection after freeing of large objects. ")
# argv = parse_args(argp, argv=c("--input", "../../../gridss-purple-linx/test/gridss/COLO829v001R_COLO829v001T.gridss.vcf", "--output", "../../../temp/somatic.vcf", "-f", "../../../temp/full.vcf", "-p", "../../../gridss-purple-linx/refdata/hg19/dbs/gridss/pon3792v1", "--scriptdir", "../", "--gc"))
# argv = parse_args(argp, argv=c("--input", "C:/dev/colo829hg38/out.vcf", "--output", "C:/temp/tmp.vcf", "-f", "C:/temp/tmp-full.vcf", "-p", "S:/hartwig/pon", "--scriptdir", "C:/dev/gridss/scripts", "--gc", "--tumourordinal", "2"))
# argv = parse_args(argp, argv=c("--input", "C:/dev/gridss/scripts/gridss2_manuscript/publicdata/sim/gen/run1_/gridss/run1_.gridss.vcf", "--output", "C:/temp/tmp.vcf", "--scriptdir", "C:/dev/gridss/scripts", "--gc", "--tumourordinal", "2"))
argv = parse_args(argp)

for (argname in c("input")) {
	if (is.na(argv[argname]) || is.null(argv[argname])) {
		msg = paste0("Required argument missing: --", argname)
		write(msg, stderr())
		print(argp)
		stop(msg)
	}
}
if (is.na(argv$input) || is.null(argv$input) || !file.exists(argv$input)) {
  msg = paste(argv$input, "not found")
  write(msg, stderr())
  print(argp)
  stop(msg)
}
if ((is.na(argv$output) || is.null(argv$output)) && is.na(argv$fulloutput) || is.null(argv$fulloutput)) {
  msg = "Must specify at least one of --output and --fulloutput"
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
if (!is.na(argv$configdir) && !is.null(argv$configdir)) {
	options(gridss.config.dir=argv$configdir)
}
options(tidyverse.quiet = TRUE)
library(tidyverse, warn.conflicts=FALSE, quietly=TRUE)
library(readr, warn.conflicts=FALSE, quietly=TRUE)
library(stringr, warn.conflicts=FALSE, quietly=TRUE)
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
refgenome = NULL
if (!is.null(argv$ref) & !is.na(argv$ref) & argv$ref != "") {
  if (!(argv$ref %in% installed.packages()[,1])) {
    stop(paste("Missing reference genome package", argv$ref, "."))
  } else {
    refgenome=eval(parse(text=paste0("library(", argv$ref, ")\n", argv$ref)))
  }
} else {
  msg = paste("No reference genome supplied using --ref. Not performing variant equivalence checks.")
  write(msg, stderr())
}
library(ggplot2)
dpi=300
theme_set(theme_bw())

# Filter to somatic calls
write(paste(Sys.time(), "Reading", argv$input), stderr())
raw_vcf = readVcf(argv$input)
nsamples = ncol(geno(raw_vcf)$VF)
if (is.null(argv$tumourordinal) | any(is.na(argv$tumourordinal))) {
	argv$tumourordinal = seq(ncol(geno(raw_vcf)$VF))[-argv$normalordinal]
}
if (any(argv$normalordinal > nsamples)) {
	msg = paste("Unable to use sample(s)", paste(argv$normalordinal, collapse = ","), "as matched normal - only", nsamples, "samples found in VCF.")
	write(msg, stderr())
	stop(msg)
}
if (any(argv$tumourordinal > nsamples)) {
	msg = paste("Unable to use sample(s)", paste(argv$tumourordinal, collapse = ","), "as tumour sample - only", nsamples, "samples found in VCF.")
	write(msg, stderr())
	stop(msg)
}

write(paste("Tumour samples:", paste(colnames(geno(raw_vcf)$VF)[argv$tumourordinal], collapse = ",")), stderr())
write(paste("Matched normals:", paste(colnames(geno(raw_vcf)$VF)[argv$normalordinal], collapse = ",")), stderr())
if (argv$plotdir != "") {
	dir.create(argv$plotdir, showWarnings=FALSE, recursive=TRUE)
  # QUAL distributions of each sample
  plotdf = data.frame(
  		sample=rep(colnames(geno(raw_vcf)$QUAL), each=length(raw_vcf)),
  		bpqual=as.numeric(geno(raw_vcf)$QUAL),
  		bequal=as.numeric(geno(raw_vcf)$BQ),
  		alt=rep(as.character(alt(raw_vcf)), times=ncol(geno(raw_vcf)$QUAL)),
  		stringsAsFactors=FALSE) %>%
		mutate(
			isSingleBreakend=startsWith(alt, ".") | endsWith(alt, "."),
			qual=ifelse(isSingleBreakend, bequal, bpqual),
			breakpointType=ifelse(isSingleBreakend, "Breakend", "Breakpoint"),
			foundIn=ifelse(isSingleBreakend,
										 rep(rowSums(geno(raw_vcf)$BQ > 0), times=ncol(geno(raw_vcf)$BQ)),
										 rep(rowSums(geno(raw_vcf)$QUAL > 0), times=ncol(geno(raw_vcf)$QUAL))))
  ggsave(paste0(argv$plotdir, "/", basename(argv$output), ".plots.bpqual.png"),
         ggplot(plotdf %>% filter(qual > 0)) +
           aes(x=qual, fill=as.character(foundIn)) +
           geom_histogram(bins=50) +
           scale_x_log10() +
           facet_wrap(sample ~ breakpointType, scales="free_y") +
           labs(fill="Found in samples", title="Raw QUAL breakdown"),
         dpi=dpi, unit="in", width=3840/dpi, height=2160/dpi)
}
# hard filter variants that are obviously not somatic
full_vcf = raw_vcf[geno(raw_vcf)$QUAL[,argv$normalordinal] / VariantAnnotation::fixed(raw_vcf)$QUAL < 4 * gridss.allowable_normal_contamination]
rm(raw_vcf)
if (argv$gc) { gc() }
# hard filter unpaired breakpoints (caused by inconsistent scoring across the two breakends)
full_vcf = full_vcf[is.na(info(full_vcf)$MATEID) | info(full_vcf)$MATEID %in% names(full_vcf)]
full_vcf = align_breakpoints(full_vcf)
# Add header fields
full_vcf = addVCFHeaders(full_vcf)

info(full_vcf)$TAF = rep("", length(full_vcf))
filters = rep("", length(full_vcf))
names(filters) = names(full_vcf)

write(paste(Sys.time(), "Parsing single breakends", argv$input), stderr())
begr = breakpointRanges(full_vcf, unpartneredBreakends=TRUE)
if (is.null(begr$sourceId) & !is.null(begr$vcfId)) {
  stop("StructuralVariantAnnotation version mismatch.")
}
write(paste(Sys.time(), "Calculating single breakend VAF", argv$input), stderr())
begr$af = round(gridss_be_af(begr, full_vcf, argv$tumourordinal), 5)
begr$af_str = as.character(begr$af)
if (length(begr) > 0) {
  info(full_vcf[names(begr)])$TAF = begr$af_str
}
write(paste(Sys.time(), "Filtering single breakends", argv$input), stderr())
befiltered = gridss_breakend_filter(begr, full_vcf, pon_dir=argv$pondir, normalOrdinal=argv$normalordinal, tumourOrdinal=argv$tumourordinal)
filters[names(begr)] = befiltered
rm(befiltered)
full_vcf = full_vcf[passes_very_hard_filters(filters)]
filters = filters[passes_very_hard_filters(filters)]
begr = begr[names(begr) %in% names(full_vcf)]
if (argv$gc) { gc() }


write(paste(Sys.time(), "Parsing breakpoints", argv$input), stderr())
bpgr = breakpointRanges(full_vcf, unpartneredBreakends=FALSE)
if (is.null(bpgr$sourceId) & !is.null(bpgr$vcfId)) {
  stop("StructuralVariantAnnotation version mismatch.")
}
write(paste(Sys.time(), "Calculating breakpoint VAF", argv$input), stderr())
bpgr$af = round(gridss_bp_af(bpgr, full_vcf, argv$tumourordinal), 5)
bpgr$af_str = paste(bpgr$af, partner(bpgr)$af, sep=",")
if (length(bpgr) > 0) {
  info(full_vcf[names(bpgr)])$TAF = bpgr$af_str
}
write(paste(Sys.time(), "Filtering breakpoints", argv$input), stderr())
bpfiltered = gridss_breakpoint_filter(bpgr, full_vcf, bsgenome=refgenome, pon_dir=argv$pondir, normalOrdinal=argv$normalordinal, tumourOrdinal=argv$tumourordinal)
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
unpaired_breakpoint = !is.na(info(full_vcf)$MATEID) & !(info(full_vcf)$MATEID %in% names(full_vcf))
full_vcf = full_vcf[!unpaired_breakpoint]
filters = filters[passes_very_hard_filters(filters)]
filters = filters[!unpaired_breakpoint]
rm(unpaired_breakpoint)

vcf = full_vcf[passes_soft_filters(filters)]
vcf = vcf[is.na(info(vcf)$MATEID) | info(vcf)$MATEID %in% names(vcf)]
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
  !((!is.na(info(vcf)$MATEID) & info(vcf)$ASSR + info(vcf)$SR + info(vcf)$IC > 0) |
      (is.na(info(vcf)$MATEID) & info(vcf)$BASSR + info(vcf)$BSC > 0))
filters[names(vcf)[is_imprecise]] = paste0(filters[names(vcf)[is_imprecise]], ";imprecise")
filters[transitive_df$linked_by] = paste0(filters[transitive_df$linked_by], ";transitive")

# write out filtering overview plots
if (argv$plotdir != "") {
  fdf = data.frame(
    id=names(filters),
    isSingleBreakend=str_detect(unlist(alt(full_vcf)), stringr::fixed(".")),
    filter=filters,
    gridssFilters=rowRanges(full_vcf)$FILTER,
    qual=rowRanges(full_vcf)$QUAL)
  fldf = fdf %>%
    mutate(nfilters=str_count(filters, ";")) %>%
    mutate(filter=str_replace(filter, "^;", "")) %>%
    separate_rows(filter, sep=";")
  ggsave(paste0(argv$plotdir, "/", basename(argv$output), ".plots.bpfilters.png"),
    ggplot(fldf %>% filter(!isSingleBreakend & qual > 250)) +
      aes(x=qual, fill=as.character(nfilters)) +
      geom_histogram(bins=50) +
      scale_x_log10() +
      facet_wrap( ~ filter) +
      labs(fill="Filters", title="Filtering breakdown (breakpoint variants)"),
    dpi=dpi, unit="in", width=3840/dpi, height=2160/dpi)
  ggsave(paste0(argv$plotdir, "/", basename(argv$output), ".plots.befilters.png"),
    ggplot(fldf %>% filter(isSingleBreakend & qual > 500)) +
      aes(x=qual, fill=as.character(nfilters)) +
      geom_histogram(bins=50) +
      scale_x_log10() +
      facet_wrap( ~ filter) +
      labs(fill="Filters", title="Filtering breakdown (single breakend variants)"),
    dpi=dpi, unit="in", width=3840/dpi, height=2160/dpi)
  ggsave(paste0(argv$plotdir, "/", basename(argv$output), ".plots.bpassembly.png"),
      ggplot(fdf %>% filter(filter=="")) +
      aes(x=qual, fill=str_replace(gridssFilters, "(PASS)|(;?LOW_QUAL;?)", "")) +
      geom_histogram(bins=100) +
      scale_x_log10() +
      labs("Assembly support of unfiltered variants", fill="Status"),
      dpi=dpi, unit="in", width=3840/dpi, height=2160/dpi)
}

vcf = full_vcf[passes_soft_filters(filters)]
vcf = vcf[is.na(info(vcf)$MATEID) | info(vcf)$MATEID %in% names(vcf)]
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
  mutate(pass=passes_final_filters(vcf[sourceId])) %>%
  group_by(linking_group) %>%
  mutate(pass=any(pass)) %>%
  ungroup() %>%
  filter(pass)

write(paste(Sys.time(),"Calculating bebe insertion links", argv$input), stderr())
# Insertion linkage
bebeins_link_df = linked_by_breakend_breakend_insertion_classification(begr) %>%
  group_by(linked_by) %>%
  mutate(pass=passes_final_filters(vcf[sourceId])) %>%
  mutate(pass=any(pass)) %>%
  ungroup() %>%
  filter(pass) %>%
  mutate(type="bebeins")
write(paste(Sys.time(),"Calculating bebp insertion links", argv$input), stderr())
bebpins_link_df = linked_by_breakpoint_breakend_insertion_classification(bpgr, begr) %>%
  group_by(linked_by) %>%
  mutate(pass=passes_final_filters(vcf[sourceId])) %>%
  mutate(pass=any(pass)) %>%
  ungroup() %>%
  filter(pass) %>%
  mutate(type="bebpins")
# Inversion linkage
write(paste(Sys.time(),"Calculating simple inversions", argv$input), stderr())
inv_link_df = linked_by_simple_inversion_classification(bpgr) %>%
  group_by(linked_by) %>%
  mutate(pass=passes_final_filters(vcf[sourceId])) %>%
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
  mutate(pass=passes_final_filters(vcf[sourceId])) %>%
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
  dplyr::select(sourceId, linked_by) %>%
	filter(sourceId %in% row.names(vcf)) %>%
  mutate(
    QUAL=rowRanges(vcf)[sourceId]$QUAL,
    hasPolyA=str_detect(as.character(unlist(rowRanges(vcf[sourceId])$ALT)), "A{16}")) %>%
  group_by(linked_by) %>%
  # filter events where supporting fragment counts differ by too much
  mutate(
    max_supporting_fragment_count = max(ifelse(is.na(info(vcf[sourceId])$MATEID), info(vcf[sourceId])$BVF, info(vcf[sourceId])$VF)),
    min_supporting_fragment_count = min(ifelse(is.na(info(vcf[sourceId])$MATEID), info(vcf[sourceId])$BVF, info(vcf[sourceId])$VF)),
    hasPolyA=any(hasPolyA)
    ) %>%
  filter(min_supporting_fragment_count >= gridss.min_rescue_portion * max_supporting_fragment_count | hasPolyA)

write(paste(Sys.time(),"Calculating final linkage annotation", argv$input), stderr())
# Only keep the best QUAL event linkage
event_link_df = event_link_df %>%
  group_by(linked_by) %>%
  mutate(linkQUAL = pmin(QUAL)) %>%
  group_by(sourceId) %>%
  filter(QUAL == linkQUAL) %>%
  group_by(linked_by)
# Don't event link to PON filtered variants
event_link_df = event_link_df %>%
  filter(!str_detect(filters[sourceId], "PON"))
# Fix up pairing
event_link_df = event_link_df %>%
  filter(n() == 2) %>%
  ungroup()

# include both breakends of any linked breakpoints
# as linkage can be breakend specific (e.g. assembly, bpbeins)
link_rescue = bind_rows(link_df, event_link_df) %>% pull(sourceId) %>% unique()
link_rescue = c(link_rescue, bpgr[link_rescue[link_rescue %in% names(bpgr)]]$partner)

# Note that we don't rescue equivalent events
begr$partner = rep(NA, length(begr))
if (!is.null(refgenome)) {
  eqv_link_df = linked_by_equivalent_variants(full_vcf, as(rbind(as.data.frame(bpgr), as.data.frame(begr)), "GRanges"), bsgenome=refgenome) %>%
    filter(passes_final_filters(vcf[sourceId]) | sourceId %in% link_rescue) %>%
    group_by(linked_by) %>%
    filter(n() == 2) %>%
    ungroup() %>%
    mutate(type="eqv")
} else {
  eqv_link_df = NULL
}

link_summary_df = bind_rows(link_df, event_link_df, eqv_link_df) %>%
  group_by(sourceId) %>%
  summarise(linked_by=paste0(linked_by, collapse=","))

# Add linking information
info(full_vcf)$LOCAL_LINKED_BY = rep("", length(full_vcf))
info(full_vcf)$REMOTE_LINKED_BY = rep("", length(full_vcf))
info(full_vcf[link_summary_df$sourceId])$LOCAL_LINKED_BY = link_summary_df$linked_by
info(full_vcf[!is.na(info(full_vcf)$MATEID)])$REMOTE_LINKED_BY = info(full_vcf[info(full_vcf[!is.na(info(full_vcf)$MATEID)])$MATEID])$LOCAL_LINKED_BY

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
  vcf = vcf[is.na(info(vcf)$MATEID) | info(vcf)$MATEID %in% names(vcf)]
  writeVcf(vcf, argv$output, index=TRUE)
}
if (!is.na(argv$fulloutput)) {
  write(paste(Sys.time(),"Writing ", argv$fulloutput), stderr())
  vcf = full_vcf[passes_very_hard_filters(filters)]
  vcf = vcf[is.na(info(vcf)$MATEID) | info(vcf)$MATEID %in% names(vcf)]
  writeVcf(vcf, argv$fulloutput, index=TRUE)
}





