package au.edu.wehi.idsv;

import com.google.common.collect.Ordering;

public abstract class DirectedEvidenceOrder {
	public static Ordering<DirectedEvidence> ByStartEnd = new Ordering<DirectedEvidence>() {
		public int compare(DirectedEvidence arg1, DirectedEvidence arg2) {
			return BreakendSummary.ByStartEnd.compare(arg1.getBreakendSummary(), arg2.getBreakendSummary());
		  }
	};
}
