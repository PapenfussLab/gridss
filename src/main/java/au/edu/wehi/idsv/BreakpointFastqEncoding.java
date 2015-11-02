package au.edu.wehi.idsv;

import com.google.common.base.Charsets;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqRecord;

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
	private static int getStartPosition(SAMRecordAssemblyEvidence evidence) {
		return getStartPosition(evidence.getSAMRecord());
	}
	private static int getStartPosition(VariantContextDirectedEvidence evidence) {
		return evidence.getStart();
	}
	public static int getStartPosition(DirectedEvidence evidence) {
		// custom start positions are required so our order matches the sorted source SAM/BAM and VCF order 
		if (evidence instanceof SoftClipEvidence) return getStartPosition((SoftClipEvidence)evidence);
		if (evidence instanceof SAMRecordAssemblyEvidence) return getStartPosition((SAMRecordAssemblyEvidence)evidence);
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
	public static int getEncodedBreakendOffset(String fastqid) {
		String[] split =  fastqid.split("#");
		if (split.length > 2) return Integer.parseInt(split[2]);
		return 0;
	}
	public static String getEncodedID(String fastqid) {
		String[] split = fastqid.split("#");
		return fastqid.substring(split[0].length() + split[1].length() + split[2].length() + 3);
	}
	public static FastqRecord getRealignmentFastq(DirectedEvidence bp) {
		byte[] sequence = bp.getBreakendSequence();
		FastqRecord fq = null;
		if (sequence != null) {
			fq = getRealignmentFastq(
					sequence,
					bp.getBreakendQuality(),
					BreakpointFastqEncoding.getReferenceIndex(bp),
					BreakpointFastqEncoding.getStartPosition(bp),
					0,
					BreakpointFastqEncoding.getID(bp));
		}
		return fq;
	}
	/**
	 * Gets the FASTQ record for realigning the given breakend sequence
	 * @param seq (partial) breakend sequence
	 * @param qual base qualities
	 * @param referenceIndex breakend contig
	 * @param startPosition breakend anchor position
	 * @param breakendOffset start offset relative to start of full breakend sequence
	 * @param evidenceID evidenceID of breakend
	 * @return FASTQ record for breakend
	 */
	public static FastqRecord getRealignmentFastq(byte[] seq, byte[] qual, int referenceIndex, int startPosition, int breakendOffset, String evidenceID) {
		String seqHeaderPrefix = getEncodedFastqID(referenceIndex, startPosition, breakendOffset, evidenceID);
		FastqRecord fq = new FastqRecord(seqHeaderPrefix, new String(seq, Charsets.US_ASCII), "", SAMUtils.phredToFastq(qual));
		return fq;
	}
	public static String getEncodedFastqID(int referenceIndex, int startPosition, int breakendOffset, String evidenceID) {
		String seqHeaderPrefix = String.format("%s#%d#%d#%s", referenceIndex, startPosition, breakendOffset, evidenceID);
		return seqHeaderPrefix;
	}
}
