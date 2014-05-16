package au.edu.wehi.socrates;

import java.util.Comparator;

/**
 * Orders evidence by ending genomic position
 *
 */
public class DirectedEvidenceEndCoordinateComparator implements Comparator<DirectedEvidence> {
	@Override
	public int compare(DirectedEvidence arg0, DirectedEvidence arg1) {
		BreakendSummary loc0 = arg0.getBreakendSummary();
		BreakendSummary loc1 = arg1.getBreakendSummary();
		int cmp = Integer.compare(loc0.referenceIndex, loc1.referenceIndex);
		if (cmp == 0) cmp = Integer.compare(loc0.end, loc1.end);
		return cmp;
	}
}
