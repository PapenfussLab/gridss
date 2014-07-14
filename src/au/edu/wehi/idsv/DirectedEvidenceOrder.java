package au.edu.wehi.idsv;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;

public abstract class DirectedEvidenceOrder {
	public static Ordering<DirectedEvidence> ByStartEnd = new Ordering<DirectedEvidence>() {
		public int compare(DirectedEvidence arg1, DirectedEvidence arg2) {
			return BreakendSummary.ByStartEnd.compare(arg1.getBreakendSummary(), arg2.getBreakendSummary());
		}
	};
	public static Ordering<DirectedEvidence> ByEndStart = new Ordering<DirectedEvidence>() {
		public int compare(DirectedEvidence arg1, DirectedEvidence arg2) {
			return BreakendSummary.ByEndStart.compare(arg1.getBreakendSummary(), arg2.getBreakendSummary());
		}
	};
	public static Ordering<DirectedEvidence> ByStartStart2EndEnd2 = new Ordering<DirectedEvidence>() {
		public int compare(DirectedEvidence arg1, DirectedEvidence arg2) {
			BreakendSummary loc1 = arg1.getBreakendSummary();
			BreakendSummary loc2 = arg2.getBreakendSummary();
			int arg1_referenceIndex2 = 0, arg2_referenceIndex2 = 0;
			int arg1_start2 = 0, arg2_start2 = 0;
			int arg1_end2 = 0, arg2_end2 = 0;
			if (loc1 instanceof BreakpointSummary) {
				BreakpointSummary bp = (BreakpointSummary)loc1;
				arg1_referenceIndex2 = bp.referenceIndex2;
				arg1_start2 = bp.start2;
				arg1_end2 = bp.end2;
			}
			if (loc2 instanceof BreakpointSummary) {
				BreakpointSummary bp = (BreakpointSummary)loc2;
				arg2_referenceIndex2 = bp.referenceIndex2;
				arg2_start2 = bp.start2;
				arg2_end2 = bp.end2;
			}
			return ComparisonChain.start()
			        .compare(loc1.referenceIndex, loc2.referenceIndex)
			        .compare(loc1.start, loc2.start)
			        .compare(loc1.end, loc2.end)
			        .compare(arg1_referenceIndex2, arg2_referenceIndex2)
			        .compare(arg1_start2, arg2_start2)
			        .compare(arg1_end2, arg2_end2)
			        .result();
		}
	};
}
