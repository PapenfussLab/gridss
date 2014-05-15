package au.edu.wehi.socrates.vcf;

import au.edu.wehi.socrates.EvidenceMetrics;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VcfConstants {
	public static final String REALIGNMENT_FAILURE = "REALNFAIL";

	//public static final String GENE_ID = "gene_id";
	public static final String TRANSCRIPT_ID = "transcript_id";
	public static final String REALIGNMENT_SCORE = "REALNNW";
	public static final String REALIGNMENT_LENGTH = "REALNLEN";
	public static final String ASSEMBLY_CONSENSUS = "CONS";
	public static final String ASSEMBLY_QUALITY = "CONSQUAL";
	public static final String ASSEMBLY_PROGRAM = "CONSASMBR";
	public static final String ASSEMBLY_CONSENSUS_READ_COUNT = "CONSRC";
	//public static final String REALIGNMENT_EXCLUDED_BASES = "remainingSoftClip";
	public static void addHeaders(VCFHeader header) {
		// Breakpoint headers
		EvidenceAttributes.addMetaDataHeaders(header);
		header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.ASSEMBLY_CONSENSUS, 1, VCFHeaderLineType.String, "Anomolous read consensus assembly sequence"));
		header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.ASSEMBLY_QUALITY, 1, VCFHeaderLineType.Float, "Anomolous read consensus assembly quality"));
		header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.ASSEMBLY_PROGRAM, 1, VCFHeaderLineType.String, "Anomolous read consensus assembly algorithm"));
		header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.ASSEMBLY_CONSENSUS_READ_COUNT, 1, VCFHeaderLineType.Integer, "Number of anomolous reads contributing to consensus assembly"));
		header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.ASSEMBLY_CONSENSUS_READ_COUNT, 1, VCFHeaderLineType.Integer, "Number of anomolous reads contributing to consensus assembly"));
		header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.REALIGNMENT_FAILURE, 1, VCFHeaderLineType.Flag, "Breakend sequence unable to be realigned to the reference"));
		
		// Retrogene headers
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.SV_TYPE);
		//header.addMetaDataLine(new VCFInfoHeaderLine(IdsvConstants.GENE_ID, 1, VCFHeaderLineType.String, "GTF gene_id of gene containing exons"));
		header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.TRANSCRIPT_ID, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "GTF transcript_id of gene containing breakpoint exons"));
		header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.REALIGNMENT_SCORE, 1, VCFHeaderLineType.Integer, "Alignment score of exon-exon junction alignment"));
		header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.REALIGNMENT_LENGTH, 1, VCFHeaderLineType.Integer, "Number of bases realigned to another exon"));
	}
}
