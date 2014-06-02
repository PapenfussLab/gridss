package au.edu.wehi.idsv;

import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqRecord;

import com.google.common.base.Charsets;

/**
 * Determines the ID encoding of a breakpoint when performing FASTQ realignment
 * @author Daniel Cameron
 *
 */
public class BreakpointFastqEncoding {
	private BreakpointFastqEncoding() { }
	public static String getFastqID(DirectedBreakpoint bp) {
		BreakendSummary loc = bp.getBreakendSummary();
		return String.format("%s#%d#%s", loc.referenceIndex, loc.start, bp.getEvidenceID());
	}
	public static int getReferenceIndex(String fastqid) {
		return Integer.parseInt(fastqid.split("#")[0]);
	}
	public static int getStartPosition(String fastqid) {
		return Integer.parseInt(fastqid.split("#")[1]);
	}
	public static String getID(String fastqid) {
		String[] split = fastqid.split("#");
		return fastqid.substring(split[0].length() + split[1].length() + 2);
	}
	public static FastqRecord getRealignmentFastq(DirectedBreakpoint bp) {
		byte[] sequence = bp.getBreakpointSequence();
		FastqRecord fq = null;
		if (sequence != null) {
			fq = new FastqRecord(
				BreakpointFastqEncoding.getFastqID(bp),
				new String(sequence, Charsets.US_ASCII),
				"",
				SAMUtils.phredToFastq(bp.getBreakpointQuality()));
		}
		return fq;
	}
}
