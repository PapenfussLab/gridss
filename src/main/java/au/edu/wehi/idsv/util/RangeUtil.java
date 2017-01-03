package au.edu.wehi.idsv.util;

import java.util.Comparator;
import java.util.Map.Entry;

import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

public class RangeUtil {
	public static RangeSet<Integer> intersect(RangeSet<Integer> rs1, RangeSet<Integer> rs2) {
		RangeSet<Integer> rs = TreeRangeSet.create();
		for (Range<Integer> r : rs1.asRanges()) {
			rs.addAll(rs2.subRangeSet(r));
		}
		return rs;
	}
	/**
	 * Adds the given object defined over the given range where the object is the 'best'
	 * (as defined by the comparator) in the given
	 * @param rs ranges to add to 
	 * @param rangeToAdd range to add object to
	 * @param objectToAdd object to add
	 * @param comparator ordering of V used to determine 'best'
	 */
	@SuppressWarnings("rawtypes")
	public static <K extends Comparable, V> void addWhereBest(RangeMap<K, V> rm, Range<K> rangeToAdd, V objectToAdd, Comparator<V> comparator) {
		RangeSet<K> definedOver = TreeRangeSet.create();
		// work out where we should insert
		definedOver.add(rangeToAdd);
		for (Entry<Range<K>, V> entry : rm.subRangeMap(rangeToAdd).asMapOfRanges().entrySet()) {
			if (comparator.compare(entry.getValue(), objectToAdd) >= 0) {
				definedOver.remove(entry.getKey());
			}
		}
		// then do the insertion
		for (Range<K> range : definedOver.asRanges()) {
			rm.put(range, objectToAdd);
		}
	}
}
