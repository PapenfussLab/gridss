package au.edu.wehi.socrates;

import net.sf.picard.fastq.FastqRecord;
import net.sf.samtools.SAMUtils;

import com.google.common.base.Charsets;

/**
 * Determines the ID encoding of a breakpoint when performing FASTQ realignment
 * @author Daniel Cameron
 *
 */
public class BreakpointFastqEncoding {
	private BreakpointFastqEncoding() { }
	public static String getFastqID(DirectedBreakpoint bp) {
		return String.format("%s#%d#%s", bp.getReferenceIndex(), bp.getWindowStart(), bp.getEvidenceID());
	}
	public static int getReferenceIndex(String fastqid) {
		return Integer.parseInt(fastqid.split("#")[0]);
	}
	public static long getPosition(String fastqid) {
		return Long.parseLong(fastqid.split("#")[1]);
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
