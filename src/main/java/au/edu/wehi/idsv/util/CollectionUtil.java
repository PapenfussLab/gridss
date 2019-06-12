package au.edu.wehi.idsv.util;

import com.google.common.base.Function;
import com.google.common.collect.Iterables;

import java.util.ArrayList;
import java.util.List;

/**
 * List helper utility functions
 * @author Daniel Cameron
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
	/**
	 * Removes the given object from the list using reference equality, not equals()
	 * @param list
	 * @param object
	 * @return
	 */
	public static <T >boolean removeByReference(ArrayList<T> list, T object) {
		if (list == null) return false;
		int size = list.size();
		for (int i = 0; i < size; i++) {
			if (list.get(i) == object) {
				list.remove(i);
				return true;
			}
		}
		return false;
	}
	/**
	 * Adds the object to the list if it does not already exist
	 * @param list
	 * @param object
	 * @return
	 */
	public static <T> boolean addUniqueByReference(ArrayList<T> list, T object) {
		if (!containsByReference(list, object)) {
			list.add(object);
			return true;
		}
		return false;
	}
	/**
	 * Adds the object to the list if it does not already exist
	 * @param list
	 * @param object
	 * @return
	 */
	public static <T> boolean containsByReference(ArrayList<T> list, T object) {
		if (list == null) return false;
		int size = list.size();
		for (int i = 0; i < size; i++) {
			if (list.get(i) == object) {
				return true;
			}
		}
		return false;
	}
}
