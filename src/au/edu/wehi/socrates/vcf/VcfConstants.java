package au.edu.wehi.socrates.vcf;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VcfConstants {
	//public static final String GENE_ID = "gene_id";
	public static final String TRANSCRIPT_ID = "transcript_id";
	public static final String REALIGNMENT_SCORE = "REALNNW";
	public static final String REALIGNMENT_LENGTH = "REALNLEN";

	//public static final String REALIGNMENT_EXCLUDED_BASES = "remainingSoftClip";
	public static void addHeaders(VCFHeader header) {
		// Metadata headers
		VcfAttributes.addHeaders(header);
		// SV headers we use
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.SV_TYPE);
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.CONFIDENCE_INTERVAL_START_POSITION);
		
		// Retrogene headers		
		//header.addMetaDataLine(new VCFInfoHeaderLine(IdsvConstants.GENE_ID, 1, VCFHeaderLineType.String, "GTF gene_id of gene containing exons"));
		header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.TRANSCRIPT_ID, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "GTF transcript_id of gene containing breakpoint exons"));
		header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.REALIGNMENT_SCORE, 1, VCFHeaderLineType.Integer, "Alignment score of exon-exon junction alignment"));
		header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.REALIGNMENT_LENGTH, 1, VCFHeaderLineType.Integer, "Number of bases realigned to another exon"));
	}
}
