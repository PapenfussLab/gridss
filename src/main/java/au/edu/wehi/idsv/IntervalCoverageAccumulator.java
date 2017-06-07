package au.edu.wehi.idsv;

import java.io.File;
import java.util.Iterator;

import au.edu.wehi.idsv.util.IntervalAccumulator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap.Entry;
import it.unimi.dsi.fastutil.objects.ObjectBidirectionalIterator;

/**
 * Calculate binned coverage
 * 
 * @author Daniel Cameron
 *
 */
public class IntervalCoverageAccumulator {
	private final CoverageCalculationMethod method;
	private final SAMSequenceDictionary dictionary;
	private final IntervalAccumulator[] coverage;
	public IntervalCoverageAccumulator(CoverageCalculationMethod method, SAMSequenceDictionary dictionary, int binWidth, Iterator<VariantContextDirectedEvidence> it) {
		this.method = method;
		this.dictionary = dictionary;
		this.coverage = initCoverage(dictionary, binWidth, it);
	}
	private static IntervalAccumulator[] initCoverage(SAMSequenceDictionary dictionary, int binWidth, Iterator<VariantContextDirectedEvidence> it) {
		IntervalAccumulator[] coverage = new IntervalAccumulator[dictionary.getSequences().size()];
		for (int i = 0; i < coverage.length; i++) {
			coverage[i] = new IntervalAccumulator(1, dictionary.getSequence(i).getSequenceLength(), binWidth);
		}
			if (it != null) {
			while (it.hasNext()) {
				VariantContextDirectedEvidence evidence = it.next();
				BreakendSummary bs = evidence.getBreakendSummary();
				// split bins at the break-end boundary
				if (bs.direction == BreakendDirection.Forward) {
					coverage[bs.referenceIndex].splitBin(bs.nominal + 1);
				} else {
					coverage[bs.referenceIndex].splitBin(bs.nominal);
				}
			}
			}
		for (int i = 0; i < coverage.length; i++) {
			coverage[i].finaliseBins();
		}
		return coverage;
	}
	public void add(SAMRecord record, ReadGcSummary summary, double readWeight) {
		switch (method) {
		case FRAGMENT:
			coverage[summary.referenceIndex].add(summary.fragmentStart, summary.fragmentEnd, readWeight);
			break;
		case READ:
			// TODO: use actual read alignment CIGAR
			coverage[summary.referenceIndex].add(record.getAlignmentStart(), record.getAlignmentEnd(), readWeight);
			break;
		}
	}
	public void writeToBed(File bed) {
		for (int i = 0; i < coverage.length; i++) {
			String refname = dictionary.getSequence(i).getSequenceName();
			ObjectBidirectionalIterator<Entry> it = coverage[i].iterator();
			while (it.hasNext()) {
				Entry e = it.next();
				int start = e.getIntKey();
				int end = start + coverage[i].getBinSize(start) - 1;
				double value = e.getDoubleValue();
				writeBed(refname, start, end, value);
			}
		}
	}
}
