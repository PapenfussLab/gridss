package au.edu.wehi.idsv.util;

import java.util.List;

import com.google.common.base.Function;
import com.google.common.collect.Iterables;

/**
 * List helper utility functions
 * @author cameron.d
 *
 */
public class CollectionUtil {
	public static <T> int sumInt(Iterable<T> list, Function<T, Integer> f) {
		return sumInt(Iterables.transform(list, f));
	}
	public static <T> int maxInt(Iterable<T> list, Function<T, Integer> f, int defaultValue) {
		return maxInt(Iterables.transform(list, f), defaultValue);
	}
	public static int maxInt(Iterable<Integer> list, int defaultValue) {
		int max = defaultValue;
		for (Integer i : list) {
			if (i != null) max = Math.max(max, i);
		}
		return max;
	}
	public static int sumInt(Iterable<Integer> list) {
		int sum = 0;
		for (Integer i : list) {
			if (i != null) sum += i;
		}
		return sum;
	}
	public static double sumDouble(List<Double> list) {
		double sum = 0;
		for (Double i : list) {
			if (i != null) sum += i;
		}
		return sum;
	}
}
