library(StructuralVariantAnnotation, quietly=TRUE)
library(VariantAnnotation, quietly=TRUE)
library(GenomicRanges, quietly=TRUE)
library(parallel, quietly = TRUE)
library(pbapply, quietly = TRUE)
library(stringr, quietly=TRUE)
library(argparser, quietly=TRUE)

options(scipen = 999)
# Create a parser object
if(!interactive()){
  argp = arg_parser("convert VCF format")
  argp = add_argument(argp, "--reference", default="", help="Reference genome to use. Must be a valid installed BSgenome package")
  argp = add_argument(argp, "--input_vcf", help="The input vcf file")
  argp = add_argument(argp, "--output_vcf", help="The output vcf file (without the .bgz suffix)")
  argp = add_argument(argp, "--n_jobs", type="integer", default=-1, help="Number of parallel jobs")
  argp = add_argument(argp, "--interval", type="character", default=NULL, help="Interval to read")
  argv = parse_args(argp)
} else {
  argv <- list(reference=reference,
               input_vcf = input_vcf, 
               output_vcf = output_vcf, 
               n_jobs = 1, 
               interval = interval)
}
# argv <- list(
#   input_vcf = "/Users/mayalevy/Downloads/gridss/NA24385_linked.vcf.bgz",
#   output_vcf = "/Users/mayalevy/Downloads/gridss/modified_vcf.vcf.gz",
#   ref = "BSgenome.Hsapiens.UCSC.hg38",
#   n_jobs = -1
# )

# Define the path to your VCF file
# gs://cromwell-backend-ultima-data-307918/cromwell-execution/SVPipeline/9aaff528-f4e7-439e-b5b5-47ee747e2515/call-GermlineLinkVariants/NA24385_linked.vcf.bgz
#vcf_file <- "/Users/mayalevy/Downloads/gridss/401882-CL10366-Z0082-CTCTGCTGTGCAATGAT_chr1_linked_orig.vcf.bgz"
#vcf_file <- "/Users/mayalevy/Downloads/gridss/diploidSV.vcf.gz"
# gsutil cp modified_vcf_wgs.vcf.bgz gs://ultimagen-users-data/maya/deepvariant/gridss/


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

if (!is.na(argv$interval)) {
  matches <- regmatches(argv$interval, regexec("^(chr[^:]+):(\\d+)-(\\d+)$", argv$interval))
  chromosome <- matches[[1]][2]
  start_pos <- as.numeric(matches[[1]][3])
  end_pos <- as.numeric(matches[[1]][4])
  gi = GRanges(chromosome, IRanges(start_pos, end_pos))
  vcf <- readVcf(argv$input_vcf,  param=gi)

} else {
  vcf <- readVcf(argv$input_vcf,  "")
}

# Read the VCF file

# Remove all SVTYPE values from the INFO field
info(vcf)$SVTYPE <- NULL

# Extract breakpoint ranges
vcf_bp <- breakpointRanges(vcf, nominalPosition=FALSE, suffix="_bp", inferMissingBreakends = FALSE)

# Determine the simple event types for each breakpoint
vcf_bp$simpleEvent <- simpleEventType(vcf_bp)

# Initialize a simpleEvent vector with NA to match the length of the original VCF
simpleEvent <- rep(NA, length(vcf))
end_positions <- rep(NA, length(vcf))
svlens <- rep(NA, length(vcf))

# Find the matching indices between the original VCF and vcf_bp
matching_indices <- match(names(vcf), names(vcf_bp))

# Assign the simpleEvent information to the matching positions
simpleEvent[!is.na(matching_indices)] <- vcf_bp$simpleEvent[na.omit(matching_indices)]
svlens[!is.na(matching_indices)] <- vcf_bp$svLen[na.omit(matching_indices)]
end_positions[!is.na(matching_indices)] <- start(vcf)[!is.na(matching_indices)] + vcf_bp$svLen[na.omit(matching_indices)]

# Create a new header line for the END field
end_header <- DataFrame(
  Number = "1",
  Type = "Integer",
  Description = "End position of the variant described in this record",
  row.names = "END"
)

# Create a new header line for the SVLEN field
svlen_header <- DataFrame(
  Number = "1",
  Type = "Integer",
  Description = "Difference in length between REF and ALT alleles",
  row.names = "SVLEN"
)

left_svinsseq_header <- DataFrame(
  Number = ".",
  Type = "String",
  Description = "Known left side of insertion for an insertion of unknown length",
  row.names = "LEFT_SVINSSEQ"
)

right_svinsseq_header <- DataFrame(
  Number = ".",
  Type = "String",
  Description = "Known rigth side of insertion for an insertion of unknown length",
  row.names = "RIGHT_SVINSSEQ"
)

# Add the new header line to the VCF header
vcf_header <- header(vcf)

info(vcf_header) <- rbind(info(vcf_header), end_header)
info(vcf_header) <- rbind(info(vcf_header), svlen_header)
info(vcf_header) <- rbind(info(vcf_header), left_svinsseq_header)
info(vcf_header) <- rbind(info(vcf_header), right_svinsseq_header)
header(vcf) <- vcf_header

# Add the simpleEvent info to the original VCF object
info(vcf)$SVTYPE <- simpleEvent
info(vcf)$END <- end_positions
info(vcf)$SVLEN <- svlens

# remove variants with svLen smaller than 30
vcf = vcf[which(is.na(svlens) | abs(svlens) >= 30)]

### Remove MATE variant - keep only one variant  only in case of DEL or INS

# Identify indices of deletions (DEL) in the VCF
del_ins_indices <- which(!is.na(info(vcf)$SVTYPE) & ((info(vcf)$SVTYPE == "DEL") | ((info(vcf)$SVTYPE == "INS"))))

# Extract MATEID for deletions
mate_ids <- info(vcf)$MATEID[del_ins_indices]

# Initialize a logical vector to keep track of variants to remove
remove_indices <- rep(FALSE, length(vcf))

# Create a set to keep track of processed MATEIDs
processed_mate_ids <- character()
count=0
# Iterate over deletion indices to remove only one of each mate pair
for (i in seq_along(del_ins_indices)) {
  current_index <- del_ins_indices[i]
  mate_id <- info(vcf)$MATEID[[current_index]]
  var_id <- names(vcf)[current_index]
  
  if (length(mate_id) > 0 && !is.na(mate_id) && var_id %in% processed_mate_ids) {
    # If the mate ID has already been processed, mark the current variant for removal
    remove_indices[current_index] <- TRUE
  } else {
    # Mark the mate ID as processed
    processed_mate_ids <- c(processed_mate_ids, mate_id)
    info(vcf)$MATEID[[current_index]]  <- NA_character_
  }
}

# Subset the VCF to remove the marked variants
vcf <- vcf[!remove_indices]

# Change SVTYPE from INV to BND for all INV variants
inv_indices <- which(!is.na(info(vcf)$SVTYPE) & info(vcf)$SVTYPE == "INV")
info(vcf)$SVTYPE[inv_indices] <- "BND"

# Collect start and end positions for DEL variants
del_indices <- which(!is.na(info(vcf)$SVTYPE) & info(vcf)$SVTYPE == "DEL")
del_lengths <- abs(info(vcf)$SVLEN[del_indices])

# Separate DEL variants based on length
short_del_indices <- del_indices[del_lengths <= 329]
short_del_lengths <- del_lengths[del_lengths <= 329]
short_del_starts <- start(vcf)[short_del_indices]
short_del_ends <- short_del_starts + short_del_lengths - 1

# Fetch sequences for short DEL variants in one call
if (length(short_del_indices) > 0){
  short_del_seqs <- getSeq(refgenome, names = seqnames(vcf)[short_del_indices], start = short_del_starts, end = short_del_ends)
} else {
  short_del_seqs <- NULL
}

# Collect indices for INS variants
ins_indices <- which(!is.na(info(vcf)$SVTYPE) & info(vcf)$SVTYPE == "INS")

# Function to process each variant
process_variant <- function(i, vcf, short_del_indices, short_del_seqs, short_del_lengths, ins_indices) {
  result <- list()
  
  if (i %in% short_del_indices) {
    del_idx <- which(short_del_indices == i)
    deletion_length <- short_del_lengths[del_idx]
    
    # Update the REF and ALT fields for short DEL variants
    full_ref_seq <- short_del_seqs[del_idx]
    result$alt <- as.character(ref(vcf)[i])
    result$ref <- DNAString(as.character(full_ref_seq))
    result$end <- start(vcf)[i] + length(result$ref)
    
  } else if (i %in% del_indices && !(i %in% short_del_indices)) {
    # Handle long DEL variants
    result$alt <- CharacterList("<DEL>")
    result$ref <- DNAString(as.character(ref(vcf)[i]))
    result$end <- start(vcf)[i] - info(vcf)$SVLEN[i]
  } else if (i %in% ins_indices) {
    # Update the REF and ALT fields for INS variants
    alt_field = as.character(alt(vcf)[i])
    if (grepl('NNNNNNNNNN', alt_field)) {
      result$alt <- CharacterList("<INS>")
      # Extract the part before 'NNNNNNNNNN'
      left_part <- sub('NNNNNNNNNN.*', '', alt_field)
      left_part <- sub('\\.$', '', left_part)
      left_part <- sub('^.*\\.', '', left_part)
      left_part <- sub('^.*\\]', '', left_part)
      result$left_svinsseq <- left_part
      # Extract the part after 'NNNNNNNNNN'
      right_part <- sub('.*NNNNNNNNNN', '', alt_field)
      right_part <- sub('^\\.', '', right_part)
      right_part <- sub('\\..*.', '', right_part)
      right_part <- sub('\\[.*.', '', right_part)
      result$right_svinsseq <- right_part
    } else {
      result$alt <- gsub("\\[chr[0-9XY]+:[0-9]+\\[|\\]chr[0-9XY]+:[0-9]+\\]", "", alt_field)
    }
    result$ref <- DNAString(as.character(ref(vcf)[i]))
    result$end <- start(vcf)[i] + length(result$ref) - 1
  } else if (is.null(info(vcf)$SVTYPE[i])) {
    print(i)
    stop(paste("I think should not reach here. this code should be deleted", i, "."))
    alt_field <- as.character(alt(vcf)[i])
    mateid_field <- info(vcf)$MATEID[i]
    if (!is.null(mateid_field) && length(mateid_field[[1]]) != 0 && (grepl(".*\\[.*\\].*", alt_field) || grepl(".*\\].*\\[.*", alt_field))) {
      # Handle breakends BND
      result$svtype <- "BND"
      alt_length <- nchar(alt_field)
      result$svlen <- alt_length
      result$end <- start(vcf)[i]
    }
    
  }
  
  return(result)
}

# Number of cores to use
if (argv$n_jobs == -1) {
  num_cores <- detectCores() - 1
} else {
  num_cores <- argv$n_jobs
}

# Process variants in parallel with progress bar using pblapply
results <- pblapply(seq_along(vcf), process_variant, vcf = vcf, short_del_indices = short_del_indices, short_del_seqs = short_del_seqs, short_del_lengths = short_del_lengths, ins_indices = ins_indices, cl = num_cores)

# Update the VCF object with results

# Assuming results from pblapply are ready

# Initialize the necessary fields
alt_updates <- as.list(alt(vcf))
ref_updates <- as.list(ref(vcf))
end_updates <- info(vcf)$END
svtype_updates <- info(vcf)$SVTYPE
svlen_updates <- info(vcf)$SVLEN
left_svinsseq_updates <- rep(NA, length(vcf))
right_svinsseq_updates <- rep(NA, length(vcf))

# Extract updates from results
for (i in seq_along(results)) {
  if (!is.null(results[[i]]$alt)) {
    alt_updates[[i]] <- results[[i]]$alt
  }
  if (!is.null(results[[i]]$ref)) {
    ref_updates[[i]] <- results[[i]]$ref
  }
  if (!is.null(results[[i]]$end)) {
    end_updates[i] <- results[[i]]$end
  }
  if (!is.null(results[[i]]$svtype)) {
    svtype_updates[i] <- results[[i]]$svtype
  }
  if (!is.null(results[[i]]$svlen)) {
    svlen_updates[i] <- results[[i]]$svlen
  }
  if (!is.null(results[[i]]$left_svinsseq)) {
    left_svinsseq_updates[i] <- results[[i]]$left_svinsseq
  }
  if (!is.null(results[[i]]$right_svinsseq)) {
    right_svinsseq_updates[i] <- results[[i]]$right_svinsseq
  }
}

ref_updates <- lapply(ref_updates, function(x) {
  if (is.null(x)) {
    return(DNAString(""))
  } else if (is(x, "DNAString") && length(x) == 0) {
    return(DNAString(""))
  } else if (is(x, "DNAStringSet")) {
    return(as(x[1], "DNAString"))
  } else if (is.character(x)) {
    return(DNAString(x))
  } else if (is.vector(x) && all(is.na(x))) {
    return(DNAString(""))
  } else {
    return(x)
  }
})


# Apply updates to the VCF object
alt(vcf) <- as(alt_updates, "CharacterList")
ref(vcf) <- DNAStringSet(unlist(ref_updates))
info(vcf)$END <- end_updates
info(vcf)$SVTYPE <- svtype_updates
info(vcf)$SVLEN <- svlen_updates
info(vcf)$LEFT_SVINSSEQ <- left_svinsseq_updates
info(vcf)$RIGHT_SVINSSEQ <- right_svinsseq_updates

## calculating (as suggested in multiple locations the approximate genotypes)
af = geno(vcf)$AF
approx_genotype <- ifelse(af < 0.1, "0/0",
                          ifelse(af < 0.8, "0/1", "1/1"))
geno(vcf)$GT <- approx_genotype
# remove variants without SVTYPE
vcf = vcf[!is.na(info(vcf)$SVTYPE)]

# Write the modified VCF file as compressed VCF
writeVcf(vcf, filename = argv$output_vcf, index = TRUE)
