package au.edu.wehi.idsv;

import java.util.stream.Stream;

import au.edu.wehi.idsv.util.IntervalUtil;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMSequenceDictionary;

public class QueryIntervalUtil {
	public static QueryInterval[] padIntervals(SAMSequenceDictionary dictionary, QueryInterval[] intervals, int padding) {
		QueryInterval[] padded = Stream.of(intervals)
				.map(qi -> new QueryInterval(qi.referenceIndex, Math.max(1, qi.start - padding), Math.min(qi.end + padding, dictionary.getSequence(qi.referenceIndex).getSequenceLength())))
				.toArray(QueryInterval[]::new);
		QueryInterval[] optimised = QueryInterval.optimizeIntervals(padded);
		return optimised;
	}
	public static boolean overlaps(QueryInterval[] intervals, int referenceIndex, int position) {
		return overlaps(intervals, referenceIndex, position, position);
	}
	public static boolean overlaps(QueryInterval[] intervals, BreakendSummary be) {
		return overlaps(intervals, be.referenceIndex, be.start, be.end);
	}
	public static boolean overlaps(QueryInterval[] intervals, int referenceIndex, int start, int end) {
		for (QueryInterval i : intervals) {
			if (referenceIndex == i.referenceIndex && IntervalUtil.overlapsClosed(start, end, i.start, i.end)) {
				return true;
			}
		}
		return false;
	}
}
