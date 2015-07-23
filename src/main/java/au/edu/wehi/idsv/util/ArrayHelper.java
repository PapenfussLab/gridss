package au.edu.wehi.idsv.util;

import java.util.Arrays;

public class ArrayHelper {
	private static int[] EMPTY_INT_ARRAY = new int[0];
	private static float[] EMPTY_FLOAT_ARRAY = new float[0];
	/**
	 * Adds the values from the second array to the first,
	 * reusing the first array if possible
	 * The result will be expanded to fit the maximum number of array elements in either array
	 */
	public static int[] add(int[] first, int[] second) {
		if (first == null && second == null) return EMPTY_INT_ARRAY;
		if (second == null) return first;
		if (first == null) return Arrays.copyOf(second, second.length);
		if (first.length < second.length) {
			first = Arrays.copyOf(first, second.length);
		}
		int nOverlap = Math.min(first.length, second.length);
		for (int i = 0; i < nOverlap; i++) {
			first[i] += second[i];
		}
		return first;
	}
	/**
	 * Adds the values from the second array to the first,
	 * reusing the first array if possible
	 * The result will be expanded to fit the maximum number of array elements in either array
	 */
	public static float[] add(float[] first, float[] second) {
		if (first == null && second == null) return EMPTY_FLOAT_ARRAY;
		if (second == null) return first;
		if (first == null) return Arrays.copyOf(second, second.length);
		if (first.length < second.length) {
			first = Arrays.copyOf(first, second.length);
		}
		int nOverlap = Math.min(first.length, second.length);
		for (int i = 0; i < nOverlap; i++) {
			first[i] += second[i];
		}
		return first;
	}
	/**
	 * Adds the given value to the given position in the array, copying and padding with zeroes
	 * is required.
	 * @param array
	 * @param position
	 * @param value
	 * @return
	 */
	public static int[] add(int[] array, int position, int value) {
		int requiredSize = position + 1;
		if (array == null) {
			array = new int[requiredSize];
		}
		if (array.length < requiredSize) {
			array = Arrays.copyOf(array, requiredSize);
		}
		array[position] += value;
		return array;
	}
	/**
	 * Adds the values from the second array to the first,
	 * reusing the first array if possible
	 * The result will be expanded to fit the maximum number of array elements in either array
	 */
	public static int[] subtract(int[] first, int[] second) {
		if (first == null && second == null) return EMPTY_INT_ARRAY;
		if (second == null) return first;
		if (first == null) first = new int[second.length];
		if (first.length < second.length) {
			first = Arrays.copyOf(first, second.length);
		}
		int nOverlap = Math.min(first.length, second.length);
		for (int i = 0; i < nOverlap; i++) {
			first[i] -= second[i];
		}
		return first;
	}
}
