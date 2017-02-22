package au.edu.wehi.idsv.vcf;

import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFSimpleHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public class VcfConstants {
	public static final String VCF42BREAKEND = ".";
	public static void addHeaders(VCFHeader header) {
		// INFO headers
		for (VcfInfoAttributes attr : VcfInfoAttributes.values()) {
			if (attr.infoHeader() != null) {
				header.addMetaDataLine(attr.infoHeader());
			}
		}
		// FORMAT headers
		for (VcfFormatAttributes attr : VcfFormatAttributes.values()) {
			if (attr.formatHeader() != null) {
				header.addMetaDataLine(attr.formatHeader());
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
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.PARTNER_BREAKEND_ID);
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.BREAKEND_EVENT_ID);
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.HOMOLOGY_LENGTH);
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.HOMOLOGY_SEQUENCE);
		header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.SOMATIC_KEY));
		
		// Simple SV headers
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.SV_LENGTH);
		header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.CONFIDENCE_INTERVAL_END_POSITION);
		header.addMetaDataLine(new VCFSimpleHeaderLine("ALT", "INV", "Inversion"));
		header.addMetaDataLine(new VCFSimpleHeaderLine("ALT", "DUP", "Duplication"));
		header.addMetaDataLine(new VCFSimpleHeaderLine("ALT", "DEL", "Deletion"));
		header.addMetaDataLine(new VCFSimpleHeaderLine("ALT", "INS", "Insertion"));
		
		// Retrogene headers		
		//header.addMetaDataLine(new VCFInfoHeaderLine(IdsvConstants.GENE_ID, 1, VCFHeaderLineType.String, "GTF gene_id of gene containing exons"));
		//header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.TRANSCRIPT_ID, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "GTF transcript_id of gene containing breakpoint exons"));
		//header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.REALIGNMENT_SCORE, 1, VCFHeaderLineType.Integer, "Alignment score of exon-exon junction alignment"));
		//header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.REALIGNMENT_LENGTH, 1, VCFHeaderLineType.Integer, "Number of bases realigned to another exon"));
	}
}
