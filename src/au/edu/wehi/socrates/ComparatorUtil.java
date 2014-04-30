package au.edu.wehi.socrates;

public class ComparatorUtil {
	/**
	 * Compares two integers
	 * @param primary0
	 * @param primary1
	 * @return ascending order Comparator result
	 */
	public static int compare(int value1, int value2) {
		if (value1 < value2) return -1;
		if (value1 == value2) return 0;
		return 1;
	}
	public static int compare(long value1, long value2) {
		if (value1 < value2) return -1;
		if (value1 == value2) return 0;
		return 1;
	}
}
