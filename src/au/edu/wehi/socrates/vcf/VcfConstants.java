package au.edu.wehi.socrates.vcf;

import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLineCount;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

public class VcfConstants {
	public static final String REFERENCE_SPANNING_READ_PAIR_COUNT = "RPC";
	public static final String REFERENCE_READ_COUNT = "RRC";
	public static final String SOFT_CLIP_READ_COUNT = "SCC";
	public static final String UNMAPPED_MATE_READ_COUNT = "OEAC";
	public static final String DISCORDANT_READ_PAIR_COUNT = "DPC";
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
		header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.REFERENCE_SPANNING_READ_PAIR_COUNT, 1, VCFHeaderLineType.Integer, "Count of reference read pairs spanning this breakpoint supporting the reference allele"));
		header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.REFERENCE_READ_COUNT, 1, VCFHeaderLineType.Integer, "Count of reference reads spanning this breakpoint supporting the reference allele"));
		header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.SOFT_CLIP_READ_COUNT, 1, VCFHeaderLineType.Integer, "Count of soft clipped reads supporting this breakpoint"));
		header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.UNMAPPED_MATE_READ_COUNT, 1, VCFHeaderLineType.Integer, "Count of reads with unmapped mates supporting this breakpoint"));
		header.addMetaDataLine(new VCFInfoHeaderLine(VcfConstants.DISCORDANT_READ_PAIR_COUNT, 1, VCFHeaderLineType.Integer, "Count of discordantly paired reads supporting this breakpoint"));
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
