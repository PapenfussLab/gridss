package au.edu.wehi.idsv.bed;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import au.edu.wehi.idsv.LinearGenomicCoordinate;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.LineIterator;

/**
 * Minimal bed wrapper retaining only interval information
 * @author Daniel Cameron
 *
 */
public class IntervalBed {
	private final SAMSequenceDictionary dictionary;
	private final LinearGenomicCoordinate linear;
	private final RangeSet<Long> intervals;
	public int size() {
		return intervals.asRanges().size();
	}
	public IntervalBed(SAMSequenceDictionary dictionary, LinearGenomicCoordinate linear, File bed) throws IOException {
		this(dictionary, linear, toRangeSet(dictionary, linear, bed));
	}
	public IntervalBed(SAMSequenceDictionary dictionary, LinearGenomicCoordinate linear) {
		this(dictionary, linear, TreeRangeSet.<Long>create());
	}
	public IntervalBed(SAMSequenceDictionary dictionary, LinearGenomicCoordinate linear, QueryInterval[] intervals) {
		this(dictionary, linear, TreeRangeSet.<Long>create());
		for (QueryInterval qi : intervals) {
			addInterval(qi);
		}
	}
	private IntervalBed(SAMSequenceDictionary dictionary, LinearGenomicCoordinate linear, RangeSet<Long> intervals) {
		this.dictionary = dictionary;
		this.linear = linear;
		this.intervals = intervals;
	}
	public static IntervalBed merge(SAMSequenceDictionary dictionary, LinearGenomicCoordinate linear, Iterable<IntervalBed> list) {
		RangeSet<Long> blacklisted = TreeRangeSet.create();
		for (IntervalBed bed : list) {
			// TODO assert dictionaries and linear coordinates match
			blacklisted.addAll(bed.intervals);
		}
		return new IntervalBed(dictionary, linear, blacklisted);
	}
	private static RangeSet<Long> toRangeSet(SAMSequenceDictionary dictionary, LinearGenomicCoordinate linear, File bed) throws IOException {
		RangeSet<Long> rs = TreeRangeSet.create();
		BEDCodec codec = new BEDCodec();
	    try (AbstractFeatureReader<BEDFeature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(bed.getAbsolutePath(), codec, false)) {
			for (BEDFeature feat : reader.iterator()) {
				String chr = feat.getContig();
				int start = feat.getStart();
				int end = feat.getEnd();
				int referenceIndex = dictionary.getSequenceIndex(chr);
				rs.add(Range.closedOpen(linear.getLinearCoordinate(referenceIndex, start), linear.getLinearCoordinate(referenceIndex, end) + 1));
			}
        }
		return rs;
	}
	public static void addInterval(SAMSequenceDictionary dictionary, LinearGenomicCoordinate linear, RangeSet<Long> blacklisted, int referenceIndex, int start, int end) {
		blacklisted.add(Range.closedOpen(linear.getLinearCoordinate(referenceIndex, start), linear.getLinearCoordinate(referenceIndex, end) + 1));
	}
	public synchronized void addInterval(int referenceIndex, int start, int end) {
		addInterval(dictionary, linear, intervals, referenceIndex, start, end);
	}
	public synchronized void addInterval(QueryInterval qi) {
		addInterval(dictionary, linear, intervals, qi.referenceIndex, qi.start, qi.end);
	}
	/**
	 * Determines whether any of the intervals overlap the given interval
	 * @param referenceIndex
	 * @param start
	 * @param end
	 * @return
	 */
	public boolean overlaps(int referenceIndex, int start, int end) {
		Range<Long> r = Range.closedOpen(linear.getLinearCoordinate(referenceIndex, start), linear.getLinearCoordinate(referenceIndex, end) + 1);
		RangeSet<Long> hits = intervals.subRangeSet(r);
		return !hits.isEmpty();
	}
	/**
	 * Removes the given set of intervals
	 * @param toRemove intervals to remove
	 */
	public void remove(IntervalBed toRemove) {
		intervals.removeAll(toRemove.intervals.asRanges());
	}
	public void write(File bed, String name) throws IOException {
		try (BufferedWriter writer = Files.newBufferedWriter(bed.toPath(), StandardCharsets.US_ASCII)) {
			writer.write(String.format("track name=\"%s\" description=\"%s\" useScore=0\n", name, name));
			for (Range<Long> r : intervals.asRanges()) {
				long lower = r.lowerEndpoint();
				long upper = r.upperEndpoint();
				int referenceIndex = linear.getReferenceIndex(lower);
				int referenceIndex2 = linear.getReferenceIndex(upper);
				assert(referenceIndex == referenceIndex2);
				int bedStart = linear.getReferencePosition(lower) - 1;
				int bedEnd = linear.getReferencePosition(upper) - 1;
				writer.write(String.format("%s\t%d\t%d\n", dictionary.getSequence(referenceIndex).getSequenceName(), bedStart, bedEnd));
			}
		}
	}
	public QueryInterval[] asQueryInterval() {
		QueryInterval[] qis = new QueryInterval[intervals.asRanges().size()];
		int i = 0;
		for (Range<Long> r : intervals.asRanges()) {
			QueryInterval qi = new QueryInterval(linear.getReferenceIndex(r.lowerEndpoint()), linear.getReferencePosition(r.lowerEndpoint()), linear.getReferencePosition(r.upperEndpoint() - 1));
			qis[i++] = qi;
			if (linear.getReferenceIndex(r.upperEndpoint() - 1) != qi.referenceIndex) {
				throw new RuntimeException("Not Yet Implemented: support for interval spaning chromosomes and unpadded LinearGenomicCoordinate lookups. This should not happen. Please raise an issue at https://github.com/PapenfussLab/gridss/issues");
			}
		}
		return qis;
	}
}
