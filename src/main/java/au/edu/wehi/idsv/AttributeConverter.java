package au.edu.wehi.idsv;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import com.google.common.primitives.UnsignedBytes;
import htsjdk.variant.vcf.VCFConstants;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.List;

/**
 * Helper functions to extract an object or list of objects of the desired type
 * Objects are expected to be a (boxed) primitive, a string, or a list or array of a single type of these
 * as per Picard SAM tag and VCF attribute persistence APIs.
 * 
 * Note that when reading numeric attributes from a VCF, Picard does not perform string conversion.
 * @author Daniel Cameron
 *
 */
public class AttributeConverter {
	/**
	 * Attempts to extract the offset-th integer from the given object
	 * @param obj
	 * @param offset
	 * @param defaultValue
	 * @return
	 */
	public static int asIntListOffset(Object obj, int offset, int defaultValue) {
		List<Integer> list = asIntList(obj);
		if (list == null || list.size() <= offset || list.get(offset) == null) return defaultValue;
		return list.get(offset);
	}
	/**
	 * Attempts to extract the offset-th double from the given object
	 * @param obj
	 * @param offset
	 * @param defaultValue
	 * @return
	 */
	public static double asDoubleListOffset(Object obj, int offset, double defaultValue) {
		List<Double> list = asDoubleList(obj);
		if (list == null || list.size() <= offset || list.get(offset) == null) return defaultValue;
		return list.get(offset);
	}
	/**
	 * Attempts to convert the given object to an integer, silently widening or narrowing as required
	 * @param obj
	 * @return
	 */
	public static Integer asInt(Object obj, int defaultValue) {
		if (obj instanceof Integer) return (Integer)obj;
		return asIntListOffset(obj, 0, defaultValue);
	}
	private static Integer asRawInt(Object obj) {
		if (obj == null) return null;
		if (obj instanceof Integer) return (Integer)obj;
		if (obj instanceof Short) return (int)(short)(Short)obj;
		if (obj instanceof Long) return (int)(long)(Long)obj;
		if (obj instanceof Byte) return UnsignedBytes.toInt((byte)(Byte)obj);
		if (obj instanceof Float) return (int)(float)(Float)obj;
		if (obj instanceof Double) return (int)(double)(Double)obj;
		return Integer.parseInt(obj.toString());
	}
	public static Double asDouble(Object obj, double defaultValue) {
		if (obj instanceof Double) return (Double)obj;
		return asDoubleListOffset(obj, 0, defaultValue);
	}
	private static Double asRawDouble(Object obj) {
		if (obj == null) return null;
		if (obj instanceof Double) return (Double)obj;
		if (obj instanceof Float) return (double)(float)(Float)obj;
		if (obj instanceof Integer) return (Double)(double)(int)(Integer)obj;
		if (obj instanceof Short) return (double)(short)(Short)obj;
		if (obj instanceof Long) return (double)(long)(Long)obj;
		if (obj instanceof Byte) return (double)UnsignedBytes.toInt((byte)(Byte)obj);
		return Double.parseDouble(obj.toString());
	}
	public static List<Integer> asIntList(Object obj) {
		return asIntList(obj, ImmutableList.<Integer>of());
	}
	/**
	 * Returns the given object as a list of integers 
	 * @param obj
	 * @param defaultValue
	 * @return
	 */
	public static List<Integer> asIntList(Object obj, List<Integer> defaultValue) {
		if (obj == null || obj.equals(VCFConstants.MISSING_VALUE_v4)) return defaultValue;
		if (obj instanceof int[]) return Ints.asList((int[])obj);
		List<Integer> result;
		if (obj.getClass().isArray()) {
			result = Lists.newArrayList();
			int len = Array.getLength(obj);
			for (int i = 0; i < len; i++) {
				result.add(asRawInt(Array.get(obj, i)));
			}
			return result;
		} else if (obj instanceof String) {
			return Lists.newArrayList(asRawInt(obj));
		} else if (obj instanceof Iterable<?>) {
			result = Lists.newArrayList();
			for (Object o : (Iterable<?>)obj) {
				result.add(asRawInt(o));
			}
			return result;
		} else {
			return Lists.newArrayList(asRawInt(obj));
		}
	}
	public static List<Double> asDoubleList(Object obj) {
		return asDoubleList(obj, ImmutableList.<Double>of());
	}
	/**
	 * Returns the given object as a list of doubles 
	 * @param obj
	 * @param defaultValue
	 * @return
	 */
	public static List<Double> asDoubleList(Object obj, List<Double> defaultValue) {
		if (obj == null || obj.equals(VCFConstants.MISSING_VALUE_v4)) return defaultValue;
		if (obj instanceof double[]) return Doubles.asList((double[])obj);
		List<Double> result;
		if (obj.getClass().isArray()) {
			result = Lists.newArrayList();
			int len = Array.getLength(obj);
			for (int i = 0; i < len; i++) {
				result.add(asRawDouble(Array.get(obj, i)));
			}
			return result;
		} else if (obj instanceof String) {
			return Lists.newArrayList(asRawDouble(obj));
		} else if (obj instanceof Iterable<?>) {
			result = Lists.newArrayList();
			for (Object o : (Iterable<?>)obj) {
				result.add(asRawDouble(o));
			}
			return result;
		} else {
			return Lists.newArrayList(asRawDouble(obj));
		}
	}
	public static List<String> asStringList(Object obj) {
		return asStringList(obj, ImmutableList.<String>of());
	}
	public static List<String> asStringList(Object obj, List<String> defaultValue) {
		if (obj == null || obj.equals(VCFConstants.MISSING_VALUE_v4)) return defaultValue;
		if (obj instanceof String) {
			return ImmutableList.of((String)obj);
		} else if (obj instanceof String[]) {
			return Arrays.asList((String[])obj);
		} else if (obj instanceof Iterable<?>) {
			List<String> list = Lists.newArrayList();
			for (Object o : (Iterable<?>)obj) {
				list.add(o.toString());
			}
			return list;
		} else {
			return ImmutableList.of(obj.toString());
		}
	}
}
