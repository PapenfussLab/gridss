package au.edu.wehi.idsv.util;

public class IntervalUtil {
	/**
	 * Determines whether the closed intervals overlap
	 * @param start1
	 * @param end1
	 * @param start2
	 * @param end2
	 * @return
	 */
	public static boolean overlapsClosed(int start1, int end1, int start2, int end2) {
		// https://fgiesen.wordpress.com/2011/10/16/checking-for-interval-overlap/
		return start1 <= end2 && start2 <= end1;
	}
	public static boolean overlapsClosed(long start1, long end1, long start2, long end2) {
		// https://fgiesen.wordpress.com/2011/10/16/checking-for-interval-overlap/
		return start1 <= end2 && start2 <= end1;
	}
	/**
	 * Width of overlap of closed intervals
	 * @param start1
	 * @param end1
	 * @param start2
	 * @param end2
	 * @return
	 */
	public static int overlapsWidthClosed(int start1, int end1, int start2, int end2) {
		if (!overlapsClosed(start1, end1, start2, end2)) return 0;
		return Math.min(end1, end2) - Math.max(start1, start2) + 1;
	}
}
