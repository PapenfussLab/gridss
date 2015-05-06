package au.edu.wehi.idsv;

import java.util.PriorityQueue;

import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.sam.SAMRecordEndCoordinateComparator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;

public class SequentialCoverageThreshold {
	private final int threshold;
	private IntervalBed bed;
	private PriorityQueue<SAMRecord> active = new PriorityQueue<SAMRecord>(new SAMRecordEndCoordinateComparator());
	private int activeIntervalReferenceIndex = -1;
	private int activeIntervalStart;
	public SequentialCoverageThreshold(SAMSequenceDictionary dictionary, LinearGenomicCoordinate linear, int thresholdCoverage) {
		if (thresholdCoverage <= 0) throw new IllegalArgumentException("Coverage threshhold must be greater than zero.");
		this.bed = new IntervalBed(dictionary, linear);
		this.threshold = thresholdCoverage;
	}
	public void acceptRecord(SAMRecord r) {
		if (r.getReadUnmappedFlag()) return;
		assert(active.peek().getReferenceIndex() <= r.getReferenceIndex());
		processBefore(r.getReferenceIndex(), r.getAlignmentStart());
		active.add(r);
		if (active.size() == threshold) {
			activeIntervalReferenceIndex = r.getReferenceIndex();
			activeIntervalStart = r.getAlignmentStart();
		}
	}
	private void processBefore(int referenceIndex, int position) {
		while (active.peek().getReferenceIndex() < referenceIndex || active.peek().getAlignmentEnd() < position) {
			if (active.size() == threshold) {
				bed.addInterval(activeIntervalReferenceIndex, activeIntervalStart, active.peek().getAlignmentEnd());
			}
			active.poll();
		}
	}
	public IntervalBed finish() {
		processBefore(Integer.MAX_VALUE, Integer.MAX_VALUE);
		return bed;
	}
}
