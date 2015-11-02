package au.edu.wehi.idsv.graph;

/**
 * Converts floating point to exact representation
 * @author Daniel Cameron
 *
 */
public class ScalingHelper {
	/**
	 * Number of discrete points per quality score.
	 * 
	 * Maximal clique calling requires exact addition and subtraction
	 * which is not possible with floating point types so weightings
	 * are converted to an exact integer type through this scaling factor.
	 */
	private static final long PER_UNIT_WEIGHT = 1 << 20;
	public static long toScaledWeight(double weight) {
		return (long)(weight * PER_UNIT_WEIGHT);
	}
	public static double toUnscaledWeight(long weight) {
		return ((double)weight) / PER_UNIT_WEIGHT;
	}
}
