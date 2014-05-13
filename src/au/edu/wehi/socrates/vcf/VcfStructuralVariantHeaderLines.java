package au.edu.wehi.socrates.vcf;

import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VcfStructuralVariantHeaderLines {
	//http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
	/** Imprecise structural variation */
	public static final VCFInfoHeaderLine IMPRECISE = new VCFInfoHeaderLine(VcfSvConstants.IMPRECISE_KEY, 0, VCFHeaderLineType.Flag, "Imprecise structural variation");
	/** Indicates a novel structural variation */
	public static final VCFInfoHeaderLine NOVEL = new VCFInfoHeaderLine(VcfSvConstants.NOVEL_KEY, 0, VCFHeaderLineType.Flag, "Indicates a novel structural variation");
	/** End position of the variant described in this record */
	public static final VCFInfoHeaderLine END = new VCFInfoHeaderLine(VcfSvConstants.END_KEY, 1, VCFHeaderLineType.Integer, "End position of the variant described in this record");
	/** Type of structural variant */
	public static final VCFInfoHeaderLine SV_TYPE = new VCFInfoHeaderLine(VcfSvConstants.SV_TYPE_KEY, 1, VCFHeaderLineType.String, "Type of structural variant");
	/** Difference in length between REF and ALT alleles */
	public static final VCFInfoHeaderLine SV_LENGTH = new VCFInfoHeaderLine(VcfSvConstants.SV_LENGTH_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Difference in length between REF and ALT alleles");
	/** Confidence interval around POS for imprecise variants */
	public static final VCFInfoHeaderLine CONFIDENCE_INTERVAL_START_POSITION = new VCFInfoHeaderLine(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY, 2, VCFHeaderLineType.Integer, "Confidence interval around POS for imprecise variants");
	/** Confidence interval around END for imprecise variants */
	public static final VCFInfoHeaderLine CONFIDENCE_INTERVAL_END_POSITION = new VCFInfoHeaderLine(VcfSvConstants.CONFIDENCE_INTERVAL_END_POSITION_KEY, 2, VCFHeaderLineType.Integer, "Confidence interval around END for imprecise variants");
	/** Confidence interval around the length of the inserted material between breakends */
	public static final VCFInfoHeaderLine CONFIDENCE_INTERVAL_LENGTH = new VCFInfoHeaderLine(VcfSvConstants.CONFIDENCE_INTERVAL_LENGTH_KEY, 2, VCFHeaderLineType.Integer, "Confidence interval around the length of the inserted material between breakends");
	/** Length of base pair identical micro-homology at event breakpoints */
	public static final VCFInfoHeaderLine HOMOLOGY_LENGTH = new VCFInfoHeaderLine(VcfSvConstants.HOMOLOGY_LENGTH_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Length of base pair identical micro-homology at event breakpoints");
	/** Sequence of base pair identical micro-homology at event breakpoints */
	public static final VCFInfoHeaderLine HOMOLOGY_SEQUENCE = new VCFInfoHeaderLine(VcfSvConstants.HOMOLOGY_SEQUENCE_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Sequence of base pair identical micro-homology at event breakpoints");
	/** ID of the assembled alternate allele in the assembly file */
	public static final VCFInfoHeaderLine BREAKPOINT_ID = new VCFInfoHeaderLine(VcfSvConstants.BREAKPOINT_ID_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "ID of the assembled alternate allele in the assembly file");
	/** Mobile element info of the form NAME,START,END,POLARITY */
	public static final VCFInfoHeaderLine MOBILE_ELEMENT_INFO = new VCFInfoHeaderLine(VcfSvConstants.MOBILE_ELEMENT_INFO_KEY, 4, VCFHeaderLineType.String, "Mobile element info of the form NAME,START,END,POLARITY");
	/** Mobile element transduction info of the form CHR,START,END,POLARITY */
	public static final VCFInfoHeaderLine MOBILE_ELEMENT_TRANSDUCTION = new VCFInfoHeaderLine(VcfSvConstants.MOBILE_ELEMENT_TRANSDUCTION_KEY, 4, VCFHeaderLineType.String, "Mobile element transduction info of the form CHR,START,END,POLARITY");
	/** ID of this element in Database of Genomic Variation */
	public static final VCFInfoHeaderLine DGVID = new VCFInfoHeaderLine(VcfSvConstants.DGVID_KEY, 1, VCFHeaderLineType.String, "ID of this element in Database of Genomic Variation");
	/** ID of this element in DBVAR */
	public static final VCFInfoHeaderLine DBVARID = new VCFInfoHeaderLine(VcfSvConstants.DBVARID_KEY, 1, VCFHeaderLineType.String, "ID of this element in DBVAR");
	/** ID of this element in DBRIP */
	public static final VCFInfoHeaderLine DBRIPID = new VCFInfoHeaderLine(VcfSvConstants.DBRIPID_KEY, 1, VCFHeaderLineType.String, "ID of this element in DBRIP");
	/** ID of mate breakends */
	public static final VCFInfoHeaderLine MATE_BREAKEND_ID = new VCFInfoHeaderLine(VcfSvConstants.MATE_BREAKEND_ID_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "ID of mate breakends");
	/** ID of partner breakend */
	public static final VCFInfoHeaderLine PARTNER_BREAKEND_ID = new VCFInfoHeaderLine(VcfSvConstants.PARTNER_BREAKEND_ID_KEY, 1, VCFHeaderLineType.String, "ID of partner breakend");
	/** ID of event associated to breakend */
	public static final VCFInfoHeaderLine BREAKEND_EVENT_ID = new VCFInfoHeaderLine(VcfSvConstants.BREAKEND_EVENT_ID_KEY, 1, VCFHeaderLineType.String, "ID of event associated to breakend");
	/** Read Depth of segment containing breakend */
	public static final VCFInfoHeaderLine BREAKEND_READ_DEPTH = new VCFInfoHeaderLine(VcfSvConstants.BREAKEND_READ_DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Read Depth of segment containing breakend");
	/** Read Depth of adjacency */
	public static final VCFInfoHeaderLine ADJACENCY_READ_DEPTH = new VCFInfoHeaderLine(VcfSvConstants.ADJACENCY_READ_DEPTH_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Read Depth of adjacency");
	/** Copy number of segment containing breakend */
	public static final VCFInfoHeaderLine BREAKEND_COPY_NUMBER = new VCFInfoHeaderLine(VcfSvConstants.BREAKEND_COPY_NUMBER_KEY, 1, VCFHeaderLineType.Integer, "Copy number of segment containing breakend");
	/** Copy number of adjacency */
	public static final VCFInfoHeaderLine ADJACENCY_COPY_NUMBER = new VCFInfoHeaderLine(VcfSvConstants.ADJACENCY_COPY_NUMBER_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Copy number of adjacency");
	/** Confidence interval around copy number for the segment */
	public static final VCFInfoHeaderLine CONFIDENCE_INTERVAL_BREAKEND_COPY_NUMBER = new VCFInfoHeaderLine(VcfSvConstants.CONFIDENCE_INTERVAL_BREAKEND_COPY_NUMBER_KEY, 2, VCFHeaderLineType.Integer, "Confidence interval around copy number for the segment");
	/** Confidence interval around copy number for the adjacency */
	public static final VCFInfoHeaderLine CONFIDENCE_INTERVAL_ADJACENCY_COPY_NUMBER = new VCFInfoHeaderLine(VcfSvConstants.CONFIDENCE_INTERVAL_ADJACENCY_COPY_NUMBER_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Confidence interval around copy number for the adjacency");
	
	/** Copy number genotype for imprecise events */
	public static final VCFFormatHeaderLine COPY_NUMBER_GENOTYPE = new VCFFormatHeaderLine(VcfSvConstants.COPY_NUMBER_GENOTYPE_KEY, 1, VCFHeaderLineType.Integer, "Copy number genotype for imprecise events");
	/** Copy number genotype quality for imprecise events */
	public static final VCFFormatHeaderLine COPY_NUMBER_GENOTYPE_QUALITY = new VCFFormatHeaderLine(VcfSvConstants.COPY_NUMBER_GENOTYPE_QUALITY_KEY, 1, VCFHeaderLineType.Float, "Copy number genotype quality for imprecise events");
	/** Copy number genotype likelihood for imprecise events */
	public static final VCFFormatHeaderLine COPY_NUMBER_GENOTYPE_LIKLIHOOD = new VCFFormatHeaderLine(VcfSvConstants.COPY_NUMBER_GENOTYPE_LIKLIHOOD_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Copy number genotype likelihood for imprecise events");
	/** Phred style probability score that the variant is novel with respect to the genome's ancestor */
	public static final VCFFormatHeaderLine NOVEL_LIKELIHOOD = new VCFFormatHeaderLine(VcfSvConstants.NOVEL_LIKELIHOOD_KEY, 1, VCFHeaderLineType.Integer, "Phred style probability score that the variant is novel with respect to the genome's ancestor");
	/** Unique haplotype identifier */
	public static final VCFFormatHeaderLine HAPLOTYPE_ID = new VCFFormatHeaderLine(VcfSvConstants.HAPLOTYPE_ID_KEY, 1, VCFHeaderLineType.Integer, "Unique haplotype identifier");
	/** Unique identifier of ancestral haplotype */
	public static final VCFFormatHeaderLine ANCESTOR_HAPLOTYPE_ID = new VCFFormatHeaderLine(VcfSvConstants.ANCESTOR_HAPLOTYPE_ID_KEY, 1, VCFHeaderLineType.Integer, "Unique identifier of ancestral haplotype");
}
