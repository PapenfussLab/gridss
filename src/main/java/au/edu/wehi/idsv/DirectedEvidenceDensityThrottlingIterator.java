package au.edu.wehi.idsv;

import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.util.DensityThrottlingIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;

/**
 * Reduces
 * @author Daniel Cameron
 *
 */
public class DirectedEvidenceDensityThrottlingIterator extends DensityThrottlingIterator<DirectedEvidence> {
	private static final Log log = Log.getInstance(DirectedEvidenceDensityThrottlingIterator.class);
	private static final boolean DEBUG_DENSITY_WIG = Boolean.valueOf(System.getProperty("gridss.writeDensityWigs", "false"));
	private final LinearGenomicCoordinate lgc;
	private final SAMSequenceDictionary dictionary;
	private final IntervalBed throttled;
	private final boolean throttleReadPairs;
	private final boolean throttleSingleReads;
	private final SAMEvidenceSource.EvidenceSortOrder iteratorSortOrder;
	private long thresholdStart = Long.MIN_VALUE;
	private int lastReferenceIndex = -1;
	private int lastPosition = -1;
	private double lastDensity = 0;
	private PrintWriter wigWriter;
	public DirectedEvidenceDensityThrottlingIterator(
			IntervalBed throttled,
			SAMSequenceDictionary dictionary,
			LinearGenomicCoordinate lgc,
			Iterator<DirectedEvidence> it,
			SAMEvidenceSource.EvidenceSortOrder iteratorSortOrder,
			int windowSize,
			double acceptDensity,
			double maxDensity,
			boolean throttleReadPairs,
			boolean throttleSingleReads){
		super(it, windowSize, acceptDensity, maxDensity);
		this.iteratorSortOrder = iteratorSortOrder;
		this.lgc = lgc;
		this.dictionary = dictionary;
		this.throttled = throttled;
		this.throttleReadPairs = throttleReadPairs;
		this.throttleSingleReads = throttleSingleReads;
	}
	@Override
	protected long getPosition(DirectedEvidence record) {
		if (iteratorSortOrder ==  SAMEvidenceSource.EvidenceSortOrder.SAMRecordStartPosition) {
			SAMRecord r = record.getUnderlyingSAMRecord();
			return lgc.getStartLinearCoordinate(r);
		} else {
			return lgc.getLinearCoordinate(record.getBreakendSummary().referenceIndex, record.getBreakendSummary().start);
		}
	}

	@Override
	protected boolean excludedFromThrottling(DirectedEvidence record) {
		return (!throttleReadPairs && record instanceof NonReferenceReadPair) ||
				(!throttleSingleReads && record instanceof SingleReadEvidence);
	}

	@Override
	public boolean hasNext() {
		boolean x = super.hasNext();
		if (!x && wigWriter != null && lastPosition != -1) {
			wigWriter.println(String.format("%d\t%.2f", lastPosition, lastDensity));
			CloserUtil.close(wigWriter);
			wigWriter = null;
		}
		return x;
	}

	@Override
	public DirectedEvidence next() {
		DirectedEvidence evidence = super.next();
		long pos = getPosition(evidence);
		if (DEBUG_DENSITY_WIG) {
			if (lgc.getReferenceIndex(pos) != lastReferenceIndex) {
				lastReferenceIndex = lgc.getReferenceIndex(pos);
				CloserUtil.close(wigWriter);
				try {
					wigWriter = new PrintWriter(File.createTempFile("assembly_density_" + dictionary.getSequence(lastReferenceIndex).getSequenceName() + "_", ".wig", new File(".")));
					wigWriter.println("variableStep	chrom=" + dictionary.getSequence(lastReferenceIndex).getSequenceName());
					lastPosition = lgc.getReferencePosition(pos);
					lastDensity =  currentDensity();
				} catch (IOException e) {
					log.error(e);
					wigWriter = null;
				}
			}
			if (lgc.getReferencePosition(pos) != lastPosition) {
				wigWriter.println(String.format("%d	%.3f", lastPosition, lastDensity));
				lastPosition = lgc.getReferencePosition(pos);
				lastDensity = currentDensity();
			}
		}
		if (!isBelowUnconditionalAcceptanceThreshold() && thresholdStart == Long.MIN_VALUE) {
			// start of new throttling region
			thresholdStart = getPosition(evidence);
		} else if (isBelowUnconditionalAcceptanceThreshold() && thresholdStart != Long.MIN_VALUE) {
			// end of throttling region
			int startReferenceIndex = lgc.getReferenceIndex(thresholdStart);
			int startPos = lgc.getReferencePosition(thresholdStart);
			int endReferenceIndex = lgc.getReferenceIndex(getPosition(evidence));
			int endPos = lgc.getReferencePosition(getPosition(evidence));
			if (startReferenceIndex == endReferenceIndex) {
				throttled.addInterval(startReferenceIndex, startPos, endPos);
			} else {
				throttled.addInterval(startReferenceIndex, startPos, dictionary.getSequence(startReferenceIndex).getSequenceLength());
				throttled.addInterval(endReferenceIndex, 1, endPos);
			}
			log.debug(String.format("Throttled assembly evidence in interval %s:%d-%d", dictionary.getSequence(startReferenceIndex).getSequenceName(), startPos, endPos));
			thresholdStart = Long.MIN_VALUE;
		}
		return evidence;
	}
}
