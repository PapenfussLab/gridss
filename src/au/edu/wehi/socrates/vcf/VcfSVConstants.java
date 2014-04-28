package au.edu.wehi.socrates.vcf;

import org.broadinstitute.variant.vcf.VCFConstants;

public class VcfSVConstants {
	/** Imprecise structural variation */
	 public static final String IMPRECISE_KEY = "IMPRECISE";
	 /** Indicates a novel structural variation */
	 public static final String NOVEL_KEY = "NOVEL";
	 /** End position of the variant described in this record */
	 public static final String END_KEY = VCFConstants.END_KEY;
	 /** Type of structural variant */
	 public static final String SV_TYPE_KEY = "SVTYPE";
	 /** Difference in length between REF and ALT alleles */
	 public static final String SV_LENGTH_KEY = "SVLEN";
	 /** Confidence interval around POS for imprecise variants */
	 public static final String CONFIDENCE_INTERVAL_START_POSITION_KEY = "CIPOS";
	 /** Confidence interval around END for imprecise variants */
	 public static final String CONFIDENCE_INTERVAL_END_POSITION_KEY = "CIEND";
	 /** Confidence interval around the length of the inserted material between breakends */
	 public static final String CONFIDENCE_INTERVAL_LENGTH_KEY = "CILEN";
	 /** Length of base pair identical micro-homology at event breakpoints */
	 public static final String HOMOLOGY_LENGTH_KEY = "HOMLEN";
	 /** Sequence of base pair identical micro-homology at event breakpoints */
	 public static final String HOMOLOGY_SEQUENCE_KEY = "HOMSEQ";
	 /** ID of the assembled alternate allele in the assembly file */
	 public static final String BREAKPOINT_ID_KEY = "BKPTID";
	 /** Mobile element info of the form NAME,START,END,POLARITY */
	 public static final String MOBILE_ELEMENT_INFO_KEY = "MEINFO";
	 /** Mobile element transduction info of the form CHR,START,END,POLARITY */
	 public static final String MOBILE_ELEMENT_TRANSDUCTION_KEY = "METRANS";
	 /** ID of this element in Database of Genomic Variation */
	 public static final String DGVID_KEY = "DGVID";
	 /** ID of this element in DBVAR */
	 public static final String DBVARID_KEY = "DBVARID";
	 /** ID of this element in DBRIP */
	 public static final String DBRIPID_KEY = "DBRIPID";
	 /** ID of mate breakends */
	 public static final String MATE_BREAKEND_ID_KEY = "MATEID";
	 /** ID of partner breakend */
	 public static final String PARTNER_BREAKEND_ID_KEY = "PARID";
	 /** ID of event associated to breakend */
	 public static final String BREAKEND_EVENT_ID_KEY = "EVENT";
	 /** Read Depth of segment containing breakend */
	 public static final String BREAKEND_READ_DEPTH_KEY = "DP";
	 /** Read Depth of adjacency */
	 public static final String ADJACENCY_READ_DEPTH_KEY = "DPADJ";
	 /** Copy number of segment containing breakend */
	 public static final String BREAKEND_COPY_NUMBER_KEY = "CN";
	 /** Copy number of adjacency */
	 public static final String ADJACENCY_COPY_NUMBER_KEY = "CNADJ";
	 /** Confidence interval around copy number for the segment */
	 public static final String CONFIDENCE_INTERVAL_BREAKEND_COPY_NUMBER_KEY = "CICN";
	 /** Confidence interval around copy number for the adjacency */
	 public static final String CONFIDENCE_INTERVAL_ADJACENCY_COPY_NUMBER_KEY = "CICNADJ";
	 /** Copy number genotype for imprecise events */
	 public static final String COPY_NUMBER_GENOTYPE_KEY = "CN";
	 /** Copy number genotype quality for imprecise events */
	 public static final String COPY_NUMBER_GENOTYPE_QUALITY_KEY = "CNQ";
	 /** Copy number genotype likelihood for imprecise events */
	 public static final String COPY_NUMBER_GENOTYPE_LIKLIHOOD_KEY = "CNL";
	 /** Phred style probability score that the variant is novel with respect to the genome's ancestor */
	 public static final String NOVEL_LIKELIHOOD_KEY = "NQ";
	 /** Unique haplotype identifier */
	 public static final String HAPLOTYPE_ID_KEY = "HAP";
	 /** Unique identifier of ancestral haplotype */
	 public static final String ANCESTOR_HAPLOTYPE_ID_KEY = "AHAP";
}
