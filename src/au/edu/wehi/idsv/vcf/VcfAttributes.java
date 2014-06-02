package au.edu.wehi.idsv.vcf;

import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.util.Arrays;

public enum VcfAttributes {
	UNMAPPED_MATE_READ_COUNT ("UMC", 1, VCFHeaderLineType.Integer, "Count of reads with unmapped mates supporting this breakpoint"),
	UNMAPPED_MATE_TOTAL_MAPQ ("UMTMQ", 1, VCFHeaderLineType.Integer, "Combined mapq of reads with unmapped mates supporting this breakpoint"),
	DISCORDANT_READ_PAIR_COUNT ("DPC", 1, VCFHeaderLineType.Integer, "Count of discordantly paired reads supporting this breakpoint"),
	DISCORDANT_READ_PAIR_TOTAL_MAPQ ("DPTMQ", 1, VCFHeaderLineType.Integer, "Combined (pair minimum) mapq of discordant read pairs supporting this breakpoint"),
	SOFT_CLIP_READ_COUNT ("SCC", 1, VCFHeaderLineType.Integer, "Count of soft clipped reads supporting this breakpoint"),
	SOFT_CLIP_TOTAL_LENGTH ("SCB", 1, VCFHeaderLineType.Integer, "Total number of soft-cliped bases supporting this breakpoint"),
	REALIGN_MAX_LENGTH ("ALNMAX", 1, VCFHeaderLineType.Integer, "Maximum number of realigned bases supporting this breakpoint"),
	REALIGN_TOTAL_LENGTH ("ALNB", 1, VCFHeaderLineType.Integer, "Total number of realigned bases supporting this breakpoint"),
	REALIGN_TOTAL_MAPQ ("ALNTMQ", 1, VCFHeaderLineType.Integer, "Total mapq of reads with realignments supporting this breakpoint"),
	ASSEMBLY_LENGTH ("ASSLEN", 1, VCFHeaderLineType.Integer, "Length of assembly supporting this breakpoint"),
	ASSEMBLY_BASES ("ASSB", 1, VCFHeaderLineType.Integer, "Total number of read bases contributing to assembly supporting this breakpoint"),
	ASSEMBLY_READS ("ASSRC", 1, VCFHeaderLineType.Integer, "Number of reads forming breakpoint assembly"),
	// --- End of summary evidence attributes *must be kept in sync with .LastEvidenceAttribute * ---
	UNMAPPED_MATE_MAX_MAPQ ("UMMQ", 1, VCFHeaderLineType.Integer, "Maximum mapq of reads with unmapped mates supporting this breakpoint"),
	DISCORDANT_READ_PAIR_MAX_MAPQ ("DPMQ", 1, VCFHeaderLineType.Integer, "Maximum (pair minimum) mapq of discordant read pairs supporting this breakpoint"),
	SOFT_CLIP_MAX_LENGTH ("SCMAX", 1, VCFHeaderLineType.Integer, "Maximum length of soft-clip supporting this breakpoint"),
	REALIGN_MAX_MAPQ ("ALNMQ", 1, VCFHeaderLineType.Integer, "Maximum mapq of reads with realignments supporting this breakpoint"),
	
	REFERENCE_SPANNING_READ_PAIR_COUNT ("RPC", 1, VCFHeaderLineType.Integer, "Count of reference read pairs spanning this breakpoint supporting the reference allele"),
	REFERENCE_READ_COUNT ("RRC", 1, VCFHeaderLineType.Integer, "Count of reference reads spanning this breakpoint supporting the reference allele"),
	
	ASSEMBLY_CONSENSUS ("ACONS", 1, VCFHeaderLineType.String, "Anomolous read consensus assembly sequence"),
	ASSEMBLY_QUALITY ("AQUAL", 1, VCFHeaderLineType.Float, "Anomolous read consensus assembly overall quality"),
	ASSEMBLY_PROGRAM ("ASMBLR", 1, VCFHeaderLineType.String, "Anomolous read consensus assembly algorithm"),
	
	REALIGNMENT_FAILURE ("REALGNFAIL", 1, VCFHeaderLineType.Flag, "Breakend sequence unable to be realigned to the reference"),
	
	CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY ("CIRPOS", 2, VCFHeaderLineType.Integer, "Confidence interval around remote breakend POS for imprecise variants");
		
	private static final VcfAttributes LastEvidenceAttribute = ASSEMBLY_READS;
    private final VCFInfoHeaderLine header;
	VcfAttributes(String name, int count, VCFHeaderLineType type, String description) {
		this.header = new VCFInfoHeaderLine(name, count, type, description);
	}
	public VCFInfoHeaderLine infoHeader() { return header; }
	public String attribute() { return header != null ? header.getID() : null; }
	
	public static boolean isEvidenceSummary(VcfAttributes attr) {
		return attr.ordinal() <= LastEvidenceAttribute.ordinal();
	}
	public static VcfAttributes[] evidenceValues() {
		return Arrays.copyOf(values(), LastEvidenceAttribute.ordinal() + 1);
	}
}
