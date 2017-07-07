package au.edu.wehi.idsv.bed;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.picard.ReferenceLookup;

/**
 * Basic data structure for a BEDPE record
 * @author Daniel Cameron
 *
 */
public class BedpeRecord {
	public final BreakpointSummary bp;
	public final String name;
	public final String[] fields;
	public double score;
	/**
	 * Parses a BEDPE record
	 * @param annotateInexactHomologyBedpe
	 * @param line
	 * @param untemplatedSequenceColumn
	 */
	public BedpeRecord(ReferenceLookup reference, String line) {
		fields = line.split("\t");
		name = fields[6];
		score = Double.parseDouble(fields[7]);
		bp = new BreakpointSummary(
				parseBreakend(reference, fields[0], fields[1], fields[2], fields[8], BreakendDirection.Forward),
				parseBreakend(reference, fields[3], fields[4], fields[5], fields[9], BreakendDirection.Backward));
	}
	private static BreakendSummary parseBreakend(ReferenceLookup reference, String chr, String strStart, String strEnd, String strDirection, BreakendDirection defaultDirection) {
		int start = Integer.parseInt(strStart) + 1; // convert back to 1-based
		int end = Integer.parseInt(strEnd);
		BreakendDirection dir = defaultDirection;
		if (strDirection.equals("+")) {
			dir = BreakendDirection.Forward;
		} else if (strDirection.equals("-")) {
			dir = BreakendDirection.Backward;
		}
		return new BreakendSummary(reference.getSequence(chr).getContigIndex(), dir, (start + end) / 2, start, end);
	}
}