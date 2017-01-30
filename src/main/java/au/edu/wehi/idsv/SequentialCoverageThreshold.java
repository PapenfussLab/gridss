package au.edu.wehi.idsv;

import au.edu.wehi.idsv.bed.IntervalBed;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import it.unimi.dsi.fastutil.longs.LongHeapPriorityQueue;

public class SequentialCoverageThreshold {
	private final int threshold;
	private final LinearGenomicCoordinate linear;
	private IntervalBed bed;
	private LongHeapPriorityQueue active = new LongHeapPriorityQueue();
	private int activeIntervalReferenceIndex = -1;
	private int activeIntervalStart;
	public SequentialCoverageThreshold(SAMSequenceDictionary dictionary, LinearGenomicCoordinate linear, int thresholdCoverage) {
		if (thresholdCoverage <= 0) throw new IllegalArgumentException("Coverage threshhold must be greater than zero.");
		this.bed = new IntervalBed(dictionary, linear);
		this.threshold = thresholdCoverage;
		this.linear = linear;
	}
	public void acceptRecord(SAMRecord r) {
		if (r.getReadUnmappedFlag()) return;
		long readStartPos = linear.getStartLinearCoordinate(r);
		long readEndPos = linear.getEndLinearCoordinate(r);
		processBefore(readStartPos);
		active.enqueue(readEndPos);
		if (active.size() == threshold) {
			activeIntervalReferenceIndex = r.getReferenceIndex();
			activeIntervalStart = r.getAlignmentStart();
		}
	}
	private void processBefore(long position) {
		while (!active.isEmpty() && active.firstLong() < position) {
			if (active.size() == threshold) {
				bed.addInterval(activeIntervalReferenceIndex, activeIntervalStart, linear.getReferencePosition(active.firstLong()));
			}
			active.dequeueLong();
		}
	}
	public IntervalBed finish() {
		processBefore(Long.MAX_VALUE);
		return bed;
	}
}
