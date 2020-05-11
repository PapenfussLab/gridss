package au.edu.wehi.idsv;

import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.util.DensityThrottlingIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;

import java.util.Iterator;

/**
 * Reduces
 * @author Daniel Cameron
 *
 */
public class SAMRecordDensityThrottlingIterator extends DensityThrottlingIterator<SAMRecord> {
	private final LinearGenomicCoordinate lgc;
	private final SAMSequenceDictionary dictionary;
	private final IntervalBed throttled;
	private SAMRecord thresholdStart = null;
	private SAMRecord thresholdEnd = null;
	public SAMRecordDensityThrottlingIterator(IntervalBed throttled, SAMSequenceDictionary dictionary, LinearGenomicCoordinate lgc, Iterator<SAMRecord> it, int windowSize, double acceptDensity, double maxDensity) {
		super(it, windowSize, acceptDensity, maxDensity);
		this.lgc = lgc;
		this.dictionary = dictionary;
		this.throttled = throttled;
	}
	@Override
	protected long getPosition(SAMRecord record) {
		return lgc.getLinearCoordinate(record.getReferenceIndex(), record.getAlignmentStart());
	}
	@Override
	protected boolean excludedFromThrottling(SAMRecord record) {
		return record.getReadUnmappedFlag();
	}

	@Override
	public boolean hasNext() {
		boolean x = super.hasNext();
		if (!x) {
			logAboveThreshold(null);
		}
		return x;
	}
	private void logAboveThreshold(SAMRecord r) {
		if (isBelowUnconditionalAcceptanceThreshold() || r == null) {
			logCurrentThreshold();
		} else {
			if (thresholdStart == null) {
				// new throttling interval
				thresholdStart = r;
			}
			thresholdEnd = r;
		}
	}
	private void logCurrentThreshold() {
		if (throttled != null && thresholdStart != null && thresholdEnd != null) {
			throttled.addInterval(thresholdStart.getReferenceIndex(), thresholdStart.getAlignmentStart(), thresholdEnd.getAlignmentStart());
		}
		thresholdStart = null;
		thresholdEnd = null;
	}

	@Override
	public SAMRecord next() {
		SAMRecord evidence = super.next();
		if (!excludedFromThrottling(evidence)) {
			logAboveThreshold(evidence);
		}
		return evidence;
	}
}
