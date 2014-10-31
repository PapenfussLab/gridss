package au.edu.wehi.idsv.vcf;

import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public class VcfConstants {
	public static final String VCF42BREAKEND = ".";
	public static final String VCF41BREAKEND_REPLACEMENT = "[<UNKNOWN>[";
	//public static final String GENE_ID = "gene_id";
	//public static final String TRANSCRIPT_ID = "transcript_id";
	//public static final String REALIGNMENT_SCORE = "REALNNW";
	//public static final String REALIGNMENT_LENGTH = "REALNLEN";
	//public static final String REALIGNMENT_EXCLUDED_BASES = "remainingSoftClip";
	public static void addHeaders(VCFHeader header) {
		// Attribute headers
		for (VcfAttributes attr : VcfAttributes.values()) {
			if (attr.infoHeader() != null) {
				header.addMetaDataLine(attr.infoHeader());
			}
		}
		// Filter headers
		for (VcfFilter filter : VcfFilter.values()) {
			if (filter.header() != null) {
				header.addMetaDataLine(filter.header());
			}
		}
		// Standard SV headers we use
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.SV_TYPE);
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.CONFIDENCE_INTERVAL_START_POSITION);
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.IMPRECISE);
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.MATE_BREAKEND_ID);
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.BREAKEND_EVENT_ID);
		header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.SOMATIC_KEY));
		
		// indel headers
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.SV_LENGTH);
		header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
		header.addMetaDataLine(new VCFInfoHeaderLine("ALT", 1, VCFHeaderLineType.String, "Temporary indel testing hack. //TODO: //FIXME: remove when testing is complete"));
		
		// Retrogene headers		
		//header.addMetaDataLine(new VCFInfoHeaderLine(IdsvConstants.GENE_ID, 1, VCFHeaderLineType.String, "GTF gene_id of gene containing exons"));
		//header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.TRANSCRIPT_ID, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "GTF transcript_id of gene containing breakpoint exons"));
		//header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.REALIGNMENT_SCORE, 1, VCFHeaderLineType.Integer, "Alignment score of exon-exon junction alignment"));
		//header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.REALIGNMENT_LENGTH, 1, VCFHeaderLineType.Integer, "Number of bases realigned to another exon"));
	}
}
