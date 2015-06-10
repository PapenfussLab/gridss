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
}
