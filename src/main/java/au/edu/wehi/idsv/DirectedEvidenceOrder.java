package au.edu.wehi.idsv;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;
import htsjdk.samtools.SAMRecord;

import java.util.Comparator;

public abstract class DirectedEvidenceOrder {
    public static final Comparator<? super DirectedEvidence> BySAMStartDeterministic = (o1, o2) -> {
		SAMRecord s1 = o1.getUnderlyingSAMRecord();
		SAMRecord s2 = o2.getUnderlyingSAMRecord();
		int cmp = new SAMRecordCoordinateOnlyComparator().fileOrderCompare(s1, s2);
		if (cmp != 0) return cmp;
		cmp = ComparisonChain.start()
				.compare(s1.getFlags(), s2.getFlags())
				.compare(s1.getReadName(), s2.getReadName())
				.result();
		if (cmp != 0) return cmp;
		// force deterministic ordering
		return o1.getEvidenceID().compareTo(o2.getEvidenceID());
	};
	public static Ordering<DirectedEvidence> ByStartEndStart2End2 = new Ordering<DirectedEvidence>() {
		public int compare(DirectedEvidence arg1, DirectedEvidence arg2) {
			BreakendSummary loc1 = arg1.getBreakendSummary();
			BreakendSummary loc2 = arg2.getBreakendSummary();
			int arg1_referenceIndex2 = 0, arg2_referenceIndex2 = 0;
			int arg1_start2 = 0, arg2_start2 = 0;
			int arg1_end2 = 0, arg2_end2 = 0;
			int arg1_nominal2 = 0, arg2_nominal2 = 0;
			if (loc1 instanceof BreakpointSummary) {
				BreakpointSummary bp = (BreakpointSummary)loc1;
				arg1_referenceIndex2 = bp.referenceIndex2;
				arg1_start2 = bp.start2;
				arg1_end2 = bp.end2;
				arg1_nominal2 = bp.nominal2;
			}
			if (loc2 instanceof BreakpointSummary) {
				BreakpointSummary bp = (BreakpointSummary)loc2;
				arg2_referenceIndex2 = bp.referenceIndex2;
				arg2_start2 = bp.start2;
				arg2_end2 = bp.end2;
				arg2_nominal2 = bp.nominal2;
			}
			int cmp = ComparisonChain.start()
			        .compare(loc1.referenceIndex, loc2.referenceIndex)
			        .compare(loc1.start, loc2.start)
			        .compare(loc1.end, loc2.end)
			        .compare(loc1.nominal, loc2.nominal)
			        .compare(arg1_referenceIndex2, arg2_referenceIndex2)
			        .compare(arg1_start2, arg2_start2)
			        .compare(arg1_end2, arg2_end2)
			        .compare(arg1_nominal2, arg2_nominal2)
			        .result();
			if (cmp != 0) return cmp;
			return BySAMStartDeterministic.compare(arg1, arg2);
		}
	};
	/**
	 * Natural (genomic location of breakend) ordering of directed evidence.  
	 */
	public static Ordering<DirectedEvidence> ByNatural = ByStartEndStart2End2;
}
