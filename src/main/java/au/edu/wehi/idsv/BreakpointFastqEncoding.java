package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqRecord;

import com.google.common.base.Charsets;

/**
 * Determines the ID encoding of a breakpoint when performing FASTQ realignment
 * 
 * @author Daniel Cameron
 *
 */
public class BreakpointFastqEncoding {
	private BreakpointFastqEncoding() { }
	public static String getID(DirectedEvidence evidence) {
		return evidence.getEvidenceID();
	}
	private static int getStartPosition(SAMRecord r) {
		return r.getAlignmentStart();
	}
	private static int getStartPosition(SoftClipEvidence evidence) {
		return getStartPosition(evidence.getSAMRecord());
	}
	private static int getStartPosition(VariantContextDirectedEvidence evidence) {
		return evidence.getStart();
	}
	public static int getStartPosition(DirectedEvidence evidence) {
		// custom start positions are required so our order matches the sorted source SAM/BAM and VCF order 
		if (evidence instanceof SoftClipEvidence) return getStartPosition((SoftClipEvidence)evidence);
		if (evidence instanceof VariantContextDirectedEvidence) return getStartPosition((VariantContextDirectedEvidence)evidence);
		return evidence.getBreakendSummary().start;
	}
	public static int getReferenceIndex(DirectedEvidence evidence) {
		return evidence.getBreakendSummary().referenceIndex;
	}
	public static int getEncodedReferenceIndex(String fastqid) {
		return Integer.parseInt(fastqid.split("#")[0]);
	}
	public static int getEncodedStartPosition(String fastqid) {
		return Integer.parseInt(fastqid.split("#")[1]);
	}
	public static String getEncodedID(String fastqid) {
		String[] split = fastqid.split("#");
		return fastqid.substring(split[0].length() + split[1].length() + 2);
	}
	public static FastqRecord getRealignmentFastq(DirectedEvidence bp) {
		byte[] sequence = bp.getBreakendSequence();
		FastqRecord fq = null;
		if (sequence != null) {
			fq = new FastqRecord(
				String.format("%s#%d#%s",
						BreakpointFastqEncoding.getReferenceIndex(bp),
						BreakpointFastqEncoding.getStartPosition(bp),
						BreakpointFastqEncoding.getID(bp)),
				new String(sequence, Charsets.US_ASCII),
				"",
				SAMUtils.phredToFastq(bp.getBreakendQuality()));
		}
		return fq;
	}
}
