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
argp = arg_parser("Link a raw GRIDSS insertions")
argp = add_argument(argp, "--ref", default="", help="Reference genome to use. Must be a valid installed BSgenome package")
argp = add_argument(argp, "--input", help="GRIDSS VCF")
argp = add_argument(argp, "--output", help="High confidence somatic subset")
argp = add_argument(argp, "--fulloutput", help="Full call set excluding obviously germline call.")
argp = add_argument(argp, "--scriptdir", default=thisFile(), help="Path to libgridss.R script")
argp = add_argument(argp, "--configdir", default=".", help="Path to gridss.config.R script relative to scriptdir. Defaults to '.' (same directory as libgridss.R)")
argp = add_argument(argp, "--gc", flag=TRUE, help="Perform garbage collection after freeing of large objects. ")
argv = parse_args(argp)

# argv <- list(
#   input = "/Users/mayalevy/Downloads/gridss/output_prefix_chr9_matt_assmebly_sv9_non_ref_trim_hap_out_sorted_ua.aln_sorted_empty.ann.vcf",
#   output = "/Users/mayalevy/Downloads/gridss/output.txt",
#   fulloutput = "/Users/mayalevy/Downloads/gridss/output_prefix_chr9_matt_assmebly_sv9_non_ref_trim_hap_out_sorted_ua.aln_sorted_empty.ann_output_of_r.vcf",
#   ref = "BSgenome.Hsapiens.UCSC.hg38",
#   scriptdir = "/Users/mayalevy/PycharmProjects/gridss/scripts/",
#   configdir = ".",
#   gc = TRUE
# )


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

raw_vcf = readVcf(argv$input)
nsamples = ncol(geno(raw_vcf)$VF)
# hard filter variants that are obviously not somatic
#full_vcf = raw_vcf[geno(raw_vcf)$QUAL[,argv$normalordinal] / VariantAnnotation::fixed(raw_vcf)$QUAL < 4 * gridss.allowable_normal_contamination]
#rm(raw_vcf)
full_vcf = raw_vcf
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
begr$af = round(gridss_be_af(begr, full_vcf, 1), 5)
begr$af_str = as.character(begr$af)
if (length(begr) > 0) {
  info(full_vcf[names(begr)])$TAF = begr$af_str
}
write(paste(Sys.time(), "Filtering single breakends", argv$input), stderr())
befiltered = gridss_breakend_filter(begr, full_vcf, pon_dir=NULL, normalOrdinal=NULL, tumourOrdinal=1, somatic_filters=FALSE)
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
bpgr$af = round(gridss_bp_af(bpgr, full_vcf, 1), 5)
bpgr$af_str = paste(bpgr$af, partner(bpgr)$af, sep=",")
if (length(bpgr) > 0) {
  info(full_vcf[names(bpgr)])$TAF = bpgr$af_str
}
write(paste(Sys.time(), "Filtering breakpoints", argv$input), stderr())
bpfiltered = gridss_breakpoint_filter(bpgr, full_vcf, bsgenome=refgenome, pon_dir=NULL, normalOrdinal=NULL, tumourOrdinal=1, somatic_filters=FALSE, min_support_filters=FALSE)
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

#vcf = full_vcf[passes_soft_filters(filters)]
vcf = full_vcf
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


#vcf = full_vcf[passes_soft_filters(filters)]
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

#remove unnecessary filters
filters <- gsub("NoAssembledRP", "", filters)

# Optionally, remove extra semicolons
filters <- gsub(";{2,}", ";", filters) # Replace multiple semicolons with a single one
filters <- gsub("^;|;$", "", filters)  # Remove leading and trailing semicolons

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
  # set MATEID be the id of the linked breakend
  vcf_ids = names(rowRanges(vcf[is.na(info(vcf)$MATEID) & (info(vcf)$LOCAL_LINKED_BY !="")]))
  vcf_linked_by = info(vcf[is.na(info(vcf)$MATEID) & (info(vcf)$LOCAL_LINKED_BY !="")])$LOCAL_LINKED_BY
  
  vcf_chr = seqnames(rowRanges(vcf[is.na(info(vcf)$MATEID) & (info(vcf)$LOCAL_LINKED_BY !="")]))
  vcf_chr = as.numeric(gsub("\\D", "", vcf_chr))
  vcf_pos = start(rowRanges(vcf[is.na(info(vcf)$MATEID) & (info(vcf)$LOCAL_LINKED_BY !="")]))
  vcf_alt_orig = rowRanges(vcf[is.na(info(vcf)$MATEID) & (info(vcf)$LOCAL_LINKED_BY !="")])$ALT
  vcf_alt = rowRanges(vcf[is.na(info(vcf)$MATEID) & (info(vcf)$LOCAL_LINKED_BY !="")])$ALT

  # Keep only the part before the comma if a comma exists
  vcf_linked_by <- str_replace(vcf_linked_by, ",.*$", "")
  mateids = rep(NA, length(vcf_ids))
  chr_pos = rep(NA, length(vcf_ids))
  
  # Order chromosomes X, Y, and M correctly if they exist
  chr_order <- c(as.character(1:22), "X", "Y", "M")
  names(chr_order) <- chr_order
  
  for (link in unique(vcf_linked_by)) {
    # Find the indices in vcf_linked_by that have this linkage value
    indices <- which(vcf_linked_by == link)
    
    # Check if exactly two indices are found
    if (length(indices) == 2) {
      # Pair the two corresponding IDs
      mateids[indices[1]] <- vcf_ids[indices[2]]
      mateids[indices[2]] <- vcf_ids[indices[1]]
      # Pair the ALT of the corresponding IDs
      # Compare chromosomes first
      first_pos = indices[1]
      second_pos = indices[2]
      if (chr_order[vcf_chr[indices[1]]] != chr_order[vcf_chr[indices[2]]]) {
        if(chr_order[vcf_chr[indices[2]]] < chr_order[vcf_chr[indices[1]]]) {
          first_pos = indices[2]
          second_pos = indices[1]
        }
      } else {
        # If chromosomes are the same, compare positions
        if (vcf_pos[indices[2]] < vcf_pos[indices[1]]) {
          first_pos = indices[2]
          second_pos = indices[1]
        }
      }
      # add the relevant part to ALT
      vcf_alt[first_pos] = paste0(vcf_alt_orig[[first_pos]], 'NNNNNNNNNN', vcf_alt_orig[[second_pos]], '[chr', vcf_chr[[second_pos]],':',vcf_pos[[second_pos]],'[')
      vcf_alt[second_pos] = paste0(']chr', vcf_chr[first_pos],':',vcf_pos[first_pos],']', vcf_alt_orig[[first_pos]], 'NNNNNNNNNN', vcf_alt_orig[[second_pos]])
    }
  }
  # add EVENT as local_linked_by
  vcf_linked_by[is.na(mateids)]=NA
  character_list <- CharacterList(vcf_alt)
  alt(vcf[is.na(info(vcf)$MATEID) & (info(vcf)$LOCAL_LINKED_BY !="")]) = character_list
  info(vcf[is.na(info(vcf)$MATEID) & (info(vcf)$LOCAL_LINKED_BY !="")])$EVENT = vcf_linked_by
  info(vcf[is.na(info(vcf)$MATEID) & (info(vcf)$LOCAL_LINKED_BY !="")])$MATEID = mateids
  
  # remove largeNoRP
  filters_before_write <- rowRanges(vcf)$FILTER
  filters_before_write <- gsub("largeNoRP;", "", filters_before_write)   # If largeNoRP is at the beginning followed by another tag.
  filters_before_write <- gsub(";largeNoRP", "", filters_before_write)   # If largeNoRP is at the end after another tag.
  filters_before_write <- gsub("^largeNoRP$", "", filters_before_write)  # If largeNoRP is the only tag.

  VariantAnnotation::fixed(vcf)$FILTER <- filters_before_write
  writeVcf(vcf, argv$fulloutput, index=TRUE)
}





