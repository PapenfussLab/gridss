package au.edu.wehi.socrates.vcf;

import java.util.HashMap;
import java.util.Map;

import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

public class VcfConstants {
	public static String VCF42BREAKEND = ".";
	public static String VCF41BREAKEND_REPLACEMENT = "[<IDSV_PLACEHOLDER>[";
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
		// Standard SV headers we use
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.SV_TYPE);
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.CONFIDENCE_INTERVAL_START_POSITION);
		
		Map<String, String> placeholderContigMap = new HashMap<String, String>();
		placeholderContigMap.put("ID", "IDSV_PLACEHOLDER"); // TODO: should this have < > ? VCF specs are unclear
		header.addMetaDataLine(new VCFContigHeaderLine(placeholderContigMap, Integer.MAX_VALUE)); // sorted last
		
		// Retrogene headers		
		//header.addMetaDataLine(new VCFInfoHeaderLine(IdsvConstants.GENE_ID, 1, VCFHeaderLineType.String, "GTF gene_id of gene containing exons"));
		//header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.TRANSCRIPT_ID, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "GTF transcript_id of gene containing breakpoint exons"));
		//header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.REALIGNMENT_SCORE, 1, VCFHeaderLineType.Integer, "Alignment score of exon-exon junction alignment"));
		//header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.REALIGNMENT_LENGTH, 1, VCFHeaderLineType.Integer, "Number of bases realigned to another exon"));
	}
}
