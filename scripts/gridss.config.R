# configuration
gridss.short_event_size_threshold = 1000

# somatic filters
gridss.allowable_normal_contamination=0.03 # matches https://github.com/hartwigmedical/pipeline/settings/include/settings.ini:BPI_CONTAMINATION_FRACTION
gridss.min_normal_depth = 8

# initial consideration filters
gridss.max_homology_length = 50
gridss.max_inexact_homology_length = 50
gridss.max_allowable_short_event_strand_bias = 0.95
gridss.single_breakend_multiplier = 1000/350 # Require this much more support for single breakends
gridss.min_af = 0.005
# How much better a single breakend call has to be for
# associated breakpoints to be considered shadow calls
gridss.shadow_breakend_multiple = 3
# final output filters
gridss.min_qual = 350

gridss.pon.min_normal_qual = 75
gridss.pon.min_samples = 1

#' Maximum gap between breakends for GRIDSS to consider the event an insertion
gridss.dsb.maxgap=35
gridss.insertion.maxgap=gridss.dsb.maxgap
gridss.inversion.maxgap=gridss.dsb.maxgap
#' Maximum size of an translocation/templated insertion
gridss.templatedinsertion.maxgap=10000

# These are excluded from the full VCF
gridss.very_hard_filters = c("normalSupport", "SRNormalSupport")
# These are included in the final VCF
gridss.soft_filters = c("PON")

# rescue limits
# fragments supporting rescued variant / fragments supporting rescuing variant
gridss.min_rescue_portion = 0.25

# consistency
gridss.min_event_size = 32
