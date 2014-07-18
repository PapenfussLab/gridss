package au.edu.wehi.idsv.vcf;

import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.util.Arrays;

public enum VcfAttributes {
	EVIDENCE_COUNT ("EC", 3, VCFHeaderLineType.Integer, "Count of  evidence supporting structural variation  (Overall,Normal,Tumour)"),
	EVIDENCE_COUNT_READPAIR ("ECRP", 3, VCFHeaderLineType.Integer, "Count of read pair evidence supporting structural variation  (Overall,Normal,Tumour)"),
	EVIDENCE_COUNT_SOFTCLIP ("ECSC", 3, VCFHeaderLineType.Integer, "Count of soft clip evidence supporting structural variation  (Overall,Normal,Tumour)"),
	EVIDENCE_COUNT_ASSEMBLY ("ECBA", 3, VCFHeaderLineType.Integer, "Count of breakend assembly evidence supporting structural variation  (Overall,Normal,Tumour)"),
	LOG_LIKELIHOOD_RATIO ("LR", 3, VCFHeaderLineType.Float, "Log-likelihood ratio of  structural variation vs reference (Overall,Normal,Tumour)"),
	LOG_LIKELIHOOD_RATIO_READPAIR ("LRRP", 3, VCFHeaderLineType.Float, "Log-likelihood ratio of read pair structural variation vs reference (Overall,Normal,Tumour)"),
	LOG_LIKELIHOOD_RATIO_SOFTCLIP ("LRSC", 3, VCFHeaderLineType.Float, "Log-likelihood ratio of soft clip structural variation vs reference (Overall,Normal,Tumour)"),
	LOG_LIKELIHOOD_RATIO_ASSEMBLY ("LRBA", 3, VCFHeaderLineType.Float, "Log-likelihood ratio of breakend assembly structural variation vs reference (Overall,Normal,Tumour)"),
	MAPPED_READPAIR ("RMRP", 3, VCFHeaderLineType.Integer, "Count of read pair evidence that maps to a remote breakend (Overall,Normal,Tumour)"),
	MAPPED_SOFTCLIP ("RMSC", 3, VCFHeaderLineType.Integer, "Count of soft clip evidence that maps to a remote breakend (Overall,Normal,Tumour)"),
	MAPPED_ASSEMBLY ("RMBA", 3, VCFHeaderLineType.Integer, "Count of breakend assembly evidence that maps to a remote breakend (Overall,Normal,Tumour)"),
	MAPQ_READPAIR_LOCAL ("MQRPL", 6, VCFHeaderLineType.Integer, "Local MAPQ of read pair evidence (Overall,Normal,Tumour, totals then maximums)"),
	MAPQ_SOFTCLIP_LOCAL ("MQSCL", 6, VCFHeaderLineType.Integer, "Local MAPQ of soft clip evidence (Overall,Normal,Tumour, totals then maximums)"),
	MAPQ_ASSEMBLY_LOCAL ("MQBAL", 6, VCFHeaderLineType.Integer, "Local MAPQ of breakend assembly evidence (Overall,Normal,Tumour, totals then maximums)"),
	MAPQ_READPAIR_REMOTE ("MQRPR", 6, VCFHeaderLineType.Integer, "Remote MAPQ of read pair evidence (Overall,Normal,Tumour, totals then maximums)"),
	MAPQ_SOFTCLIP_REMOTE ("MQSCR", 6, VCFHeaderLineType.Integer, "Remote MAPQ of soft clip evidence (Overall,Normal,Tumour, totals then maximums)"),
	MAPQ_ASSEMBLY_REMOTE ("MQBAR", 6, VCFHeaderLineType.Integer, "Remote MAPQ of breakend assembly evidence (Overall,Normal,Tumour, totals then maximums)"),
	LENGTH_READPAIR_LOCAL ("BLRPL", 6, VCFHeaderLineType.Integer, "Local length (in bases) of read pair evidence (Overall,Normal,Tumour, totals then maximums)"),
	LENGTH_SOFTCLIP_LOCAL ("BLSCL", 6, VCFHeaderLineType.Integer, "Local length (in bases) of soft clip evidence (Overall,Normal,Tumour, totals then maximums)"),
	LENGTH_ASSEMBLY_LOCAL ("BLBAL", 6, VCFHeaderLineType.Integer, "Local length (in bases) of breakend assembly evidence (Overall,Normal,Tumour, totals then maximums)"),
	LENGTH_READPAIR_REMOTE ("BLRPR", 6, VCFHeaderLineType.Integer, "Remote length (in bases) of read pair evidence (Overall,Normal,Tumour, totals then maximums)"),
	LENGTH_SOFTCLIP_REMOTE ("BLSCR", 6, VCFHeaderLineType.Integer, "Remote length (in bases) of soft clip evidence (Overall,Normal,Tumour, totals then maximums)"),
	LENGTH_ASSEMBLY_REMOTE ("BLBAR", 6, VCFHeaderLineType.Integer, "Remote length (in bases) of breakend assembly evidence (Overall,Normal,Tumour, totals then maximums)"),
	BASE_COUNT_READPAIR_LOCAL ("BCRPL", 6, VCFHeaderLineType.Integer, "Local number of bases read pair evidence (Overall,Normal,Tumour, totals then maximums)"),
	BASE_COUNT_SOFTCLIP_LOCAL ("BCSCL", 6, VCFHeaderLineType.Integer, "Local number of bases soft clip evidence (Overall,Normal,Tumour, totals then maximums)"),
	BASE_COUNT_ASSEMBLY_LOCAL ("BCBAL", 6, VCFHeaderLineType.Integer, "Local number of bases breakend assembly evidence (Overall,Normal,Tumour, totals then maximums)"),
	BASE_COUNT_READPAIR_REMOTE ("BCRPR", 6, VCFHeaderLineType.Integer, "Remote number of bases read pair evidence (Overall,Normal,Tumour, totals then maximums)"),
	BASE_COUNT_SOFTCLIP_REMOTE ("BCSCR", 6, VCFHeaderLineType.Integer, "Remote number of bases soft clip evidence (Overall,Normal,Tumour, totals then maximums)"),
	BASE_COUNT_ASSEMBLY_REMOTE ("BCBAR", 6, VCFHeaderLineType.Integer, "Remote number of bases breakend assembly evidence (Overall,Normal,Tumour, totals then maximums)"),
	BASE_QUAL_READPAIR_LOCAL ("BQRPL", 6, VCFHeaderLineType.Integer, "Local base quality of read pair evidence (Overall,Normal,Tumour, totals then maximums)"),
	BASE_QUAL_SOFTCLIP_LOCAL ("BQSCL", 6, VCFHeaderLineType.Integer, "Local base quality of soft clip evidence (Overall,Normal,Tumour, totals then maximums)"),
	BASE_QUAL_ASSEMBLY_LOCAL ("BQBAL", 6, VCFHeaderLineType.Integer, "Local base quality of breakend assembly evidence (Overall,Normal,Tumour, totals then maximums)"),
	BASE_QUAL_READPAIR_REMOTE ("BQRPR", 6, VCFHeaderLineType.Integer, "Remote base quality of read pair evidence (Overall,Normal,Tumour, totals then maximums)"),
	BASE_QUAL_SOFTCLIP_REMOTE ("BQSCR", 6, VCFHeaderLineType.Integer, "Remote base quality of soft clip evidence (Overall,Normal,Tumour, totals then maximums)"),
	BASE_QUAL_ASSEMBLY_REMOTE ("BQBAR", 6, VCFHeaderLineType.Integer, "Remote base quality of breakend assembly evidence (Overall,Normal,Tumour, totals then maximums)"),

	// --- End of summary evidence attributes *must be kept in sync with .LastEvidenceAttribute * ---
	
	REFERENCE_SPANNING_READ_PAIR_COUNT ("RPC", 3, VCFHeaderLineType.Integer, "Count of reference read pairs spanning this breakpoint supporting the reference allele (Overall,Normal,Tumour)"),
	REFERENCE_READ_COUNT ("RRC", 3, VCFHeaderLineType.Integer, "Count of reference reads spanning this breakpoint supporting the reference allele (Overall,Normal,Tumour)"),
	
	ASSEMBLY_CONSENSUS ("ACONS", 1, VCFHeaderLineType.String, "Anomolous read consensus assembly sequence"),
	ASSEMBLY_QUALITY ("AQUAL", 1, VCFHeaderLineType.Float, "Anomolous read consensus assembly overall quality"),
	ASSEMBLY_PROGRAM ("ASMBLR", 1, VCFHeaderLineType.String, "Anomolous read consensus assembly algorithm"),
	ASSEMBLY_MAX_READ_SOFT_CLIP ("ASCMAX", 1, VCFHeaderLineType.Integer, "Longest soft clip length of reads contributing to the assembly"),
	
	CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY ("CIRPOS", 2, VCFHeaderLineType.Integer, "Confidence interval around remote breakend POS for imprecise variants");
		
	private static final VcfAttributes LastEvidenceAttribute = BASE_QUAL_ASSEMBLY_REMOTE;
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
