package au.edu.wehi.socrates;

/**
 * Determines the ID encoding of a breakpoint when performing FASTQ realignment
 * @author Daniel Cameron
 *
 */
public class BreakpointFastqIDEncoding {
	private BreakpointFastqIDEncoding() { }
	public static String getFastqID(DirectedBreakpoint bp) {
		return String.format("%s#%d#%s", bp.getChr(), bp.getStart(), bp.getBreakpointID());
	}
	public static String getChr(String fastqid) {
		return fastqid.split("#")[0];
	}
	public static long getPosition(String fastqid) {
		return Long.parseLong(fastqid.split("#")[1]);
	}
	public static String getID(String fastqid) {
		String[] split = fastqid.split("#");
		return fastqid.substring(split[0].length() + split[1].length() + 2);
	}
}
