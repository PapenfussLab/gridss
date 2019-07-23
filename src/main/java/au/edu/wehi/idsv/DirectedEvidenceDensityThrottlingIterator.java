package au.edu.wehi.idsv;

import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.util.DensityThrottlingIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.FileNotFoundException;
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
	private DirectedEvidence tresholdStart = null;
	private int lastReferenceIndex = -1;
	private int lastPosition = -1;
	private double lastDensity = 0;
	private PrintWriter wigWriter;
	public DirectedEvidenceDensityThrottlingIterator(
			IntervalBed throttled,
			SAMSequenceDictionary dictionary,
			LinearGenomicCoordinate lgc,
			Iterator<DirectedEvidence> it,
			int windowSize,
			double acceptDensity,
			double maxDensity,
			boolean throttleReadPairs,
			boolean throttleSingleReads){
		super(it, windowSize, acceptDensity, maxDensity);
		this.lgc = lgc;
		this.dictionary = dictionary;
		this.throttled = throttled;
		this.throttleReadPairs = throttleReadPairs;
		this.throttleSingleReads = throttleSingleReads;
	}
	@Override
	protected long getPosition(DirectedEvidence record) {
		return lgc.getLinearCoordinate(record.getBreakendSummary().referenceIndex, record.getBreakendSummary().start);
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
		if (DEBUG_DENSITY_WIG) {
			if (evidence.getBreakendSummary().referenceIndex != lastReferenceIndex) {
				lastReferenceIndex = evidence.getBreakendSummary().referenceIndex;
				CloserUtil.close(wigWriter);
				try {
					wigWriter = new PrintWriter(File.createTempFile("assembly_density_" + dictionary.getSequence(lastReferenceIndex).getSequenceName() + "_", ".wig", new File(".")));
					wigWriter.println("variableStep	chrom=" + dictionary.getSequence(lastReferenceIndex).getSequenceName());
					lastPosition = evidence.getBreakendSummary().start;
					lastDensity =  currentDensity();
				} catch (IOException e) {
					log.error(e);
					wigWriter = null;
				}
			}
			if (evidence.getBreakendSummary().start != lastPosition) {
				wigWriter.println(String.format("%d	%.3f", lastPosition, lastDensity));
				lastPosition = evidence.getBreakendSummary().start;
				lastDensity = currentDensity();
			}
		}
		if (!isBelowUnconditionalAcceptanceThreshold() && tresholdStart == null) {
			tresholdStart = evidence;
		} else if (isBelowUnconditionalAcceptanceThreshold() && tresholdStart != null) {			
			int startReferenceIndex = tresholdStart.getBreakendSummary().referenceIndex;
			int startPos = tresholdStart.getBreakendSummary().start;
			int endReferenceIndex = evidence.getBreakendSummary().referenceIndex;
			int endPos = evidence.getBreakendSummary().start;
			if (startReferenceIndex == endReferenceIndex) {
				throttled.addInterval(startReferenceIndex, startPos, endPos);
			} else {
				throttled.addInterval(startReferenceIndex, startPos, dictionary.getSequence(startReferenceIndex).getSequenceLength());
				throttled.addInterval(endReferenceIndex, 1, endPos);
			}
			log.debug(String.format("Throttled assembly evidence in interval %s:%d-%d", dictionary.getSequence(startReferenceIndex).getSequenceName(), startPos, endPos));
			tresholdStart = null;
		}
		return evidence;
	}
}
