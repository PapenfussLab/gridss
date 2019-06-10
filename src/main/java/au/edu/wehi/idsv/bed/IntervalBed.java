package au.edu.wehi.idsv.bed;

import au.edu.wehi.idsv.LinearGenomicCoordinate;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.util.Log;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.LineIterator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.stream.Collectors;

/**
 * Minimal bed wrapper retaining only interval information
 * @author Daniel Cameron
 *
 */
public class IntervalBed {
	private static final Log log = Log.getInstance(IntervalBed.class);
	private final LinearGenomicCoordinate linear;
	private final RangeSet<Long> intervals;
	public int size() {
		return intervals.asRanges().size();
	}
	public IntervalBed(LinearGenomicCoordinate linear, File bed) throws IOException {
		this(linear, toRangeSet(linear, bed));
	}
	public IntervalBed(LinearGenomicCoordinate linear) {
		this(linear, TreeRangeSet.<Long>create());
	}
	public IntervalBed(LinearGenomicCoordinate linear, QueryInterval[] intervals) {
		this(linear, TreeRangeSet.<Long>create());
		for (QueryInterval qi : intervals) {
			addInterval(qi);
		}
	}
	private IntervalBed(LinearGenomicCoordinate linear, RangeSet<Long> intervals) {
		this.linear = linear;
		this.intervals = intervals;
	}
	public static IntervalBed merge(LinearGenomicCoordinate linear, Iterable<IntervalBed> list) {
		RangeSet<Long> blacklisted = TreeRangeSet.create();
		for (IntervalBed bed : list) {
			// TODO assert dictionaries and linear coordinates match
			blacklisted.addAll(bed.intervals);
		}
		return new IntervalBed(linear, blacklisted);
	}
	private static RangeSet<Long> toRangeSet(LinearGenomicCoordinate linear, File bed) throws IOException {
		RangeSet<Long> rs = TreeRangeSet.create();
		BEDCodec codec = new BEDCodec();
	    try (AbstractFeatureReader<BEDFeature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(bed.getPath(), codec, false)) {
	    	int lineno = 0;
			for (BEDFeature feat : reader.iterator()) {
				lineno++;
				String chr = feat.getContig();
				int start = feat.getStart();
				int end = feat.getEnd();
				int referenceIndex = linear.getDictionary().getSequenceIndex(chr);
				if (referenceIndex < 0) {
					String msg = String.format("Error loading record %d of %s: '%s' not found in reference genome.");
					log.error(msg);
					throw new IllegalArgumentException(msg);
				}
				rs.add(Range.closedOpen(linear.getLinearCoordinate(referenceIndex, start), linear.getLinearCoordinate(referenceIndex, end) + 1));
			}
        }
		return rs;
	}
	public static void addInterval(LinearGenomicCoordinate linear, RangeSet<Long> blacklisted, int referenceIndex, int start, int end) {
		blacklisted.add(Range.closedOpen(linear.getLinearCoordinate(referenceIndex, start), linear.getLinearCoordinate(referenceIndex, end) + 1));
	}
	public synchronized void addInterval(int referenceIndex, int start, int end) {
		addInterval(linear, intervals, referenceIndex, start, end);
	}
	public synchronized void addInterval(QueryInterval qi) {
		addInterval(linear, intervals, qi.referenceIndex, qi.start, qi.end);
	}
	/**
	 * Determines whether any of the intervals overlap the given interval
	 * @param referenceIndex
	 * @param start
	 * @param end
	 * @return
	 */
	public boolean overlaps(int referenceIndex, int start, int end) {
		return overlaps(linear.getLinearCoordinate(referenceIndex, start), linear.getLinearCoordinate(referenceIndex, end));
	}
	public boolean overlaps(long start, long end) {
		Range<Long> r = Range.closedOpen(start, end + 1);
		return overlaps(r);
	}
	public boolean overlaps(Range<Long> interval) {
		if (interval == null) {
			return false;
		}
		return intervals.intersects(interval);
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
				writer.write(String.format("%s\t%d\t%d\n", linear.getDictionary().getSequence(referenceIndex).getSequenceName(), bedStart, bedEnd));
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

	/**
	 * Generates a new IntervalBed in which each interval in the set by the given number of bases.
	 * Expanded intervals are truncated at reference contig bounds.
	 */
	public IntervalBed expandIntervals(int startBases, int endBases) {
		TreeRangeSet<Long> expanded = TreeRangeSet.create(intervals
				.asRanges()
				.stream()
				.map(r -> expand(r, startBases, endBases))
				.collect(Collectors.toList()));
		return new IntervalBed(linear, expanded);
	}
	private Range<Long> expand(Range<Long> range, int startBases, int endBases) {
		int referenceIndex = linear.getReferenceIndex(range.lowerEndpoint());
		int start = linear.getReferencePosition(range.lowerEndpoint());
		int end = linear.getReferencePosition(range.upperEndpoint());
		start = Math.max(1, start - startBases);
		end = Math.min(linear.getDictionary().getSequence(referenceIndex).getSequenceLength() + 1, end + endBases);
		return Range.closedOpen(linear.getLinearCoordinate(referenceIndex, start), linear.getLinearCoordinate(referenceIndex, end));
	}
	public RangeSet<Long> asRangeSet() {
		return TreeRangeSet.create(intervals);
	}
}
