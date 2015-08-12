package au.edu.wehi.idsv.sam;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

public class CigarUtil {
	private CigarUtil() {
	}
	/**
	 * Converts a xPxNxP placeholder CIGAR to the negative deletion it is intended to correspond to 
	 * @param list
	 * @return
	 */
	public static List<CigarElement> decodeNegativeDeletion(List<CigarElement> list) {
		list = Lists.newCopyOnWriteArrayList(list);
		for (int i = 0; i < list.size() - 2; i++) {
			if (list.get(i).getOperator() == CigarOperator.P &&
					list.get(i + 1).getOperator() == CigarOperator.N &&
					list.get(i + 2).getOperator() == CigarOperator.P &&
					list.get(i).getLength() == list.get(i + 1).getLength() &&
					list.get(i + 1).getLength() == list.get(i + 2).getLength()) {
				list.add(i, new CigarElement(-list.get(i).getLength(), CigarOperator.DELETION));
				list.remove(i + 1);
				list.remove(i + 1);
				list.remove(i + 1);
			}
		}
		return list;
	}
	/**
	 * Encodes a negative deletion through a xPxNxP placeholder CIGAR that conforms to the SAM specs 
	 * @param list
	 * @return valid cigar
	 */
	public static List<CigarElement> encodeNegativeDeletion(List<CigarElement> list) {
		list = Lists.newCopyOnWriteArrayList(list);
		for (int i = 0; i < list.size(); i++) {
			CigarElement e = list.get(i);
			if (e.getOperator() == CigarOperator.D && e.getLength() < 0) {
				list.remove(i);
				list.add(i, new CigarElement(-e.getLength(), CigarOperator.P));
				list.add(i, new CigarElement(-e.getLength(), CigarOperator.N));
				list.add(i, new CigarElement(-e.getLength(), CigarOperator.P));
			}
		}
		return list;
	}
	/**
	 * Returns the number of read bases 'consumed' by the given cigar 
	 * @param list
	 * @return base count
	 */
	public static int readLength(List<CigarElement> list) {
		int length = 0;
		for (CigarElement e : list) {
			if (e.getOperator().consumesReadBases()) {
				length += e.getLength();
			}
		}
		return length;
	}
	/**
	 * Returns the number of reference bases 'consumed' by the given cigar 
	 * @param list
	 * @return base count
	 */
	public static int referenceLength(List<CigarElement> list) {
		int length = 0;
		for (CigarElement e : list) {
			if (e.getOperator().consumesReferenceBases()) {
				length += e.getLength();
			}
		}
		return length;
	}
	/**
	 * Contains indel
	 * @param list
	 * @return
	 */
	public static boolean containsIndel(Cigar c) {
		for (CigarElement e : decodeNegativeDeletion(c.getCigarElements())) {
			if (e.getOperator() == CigarOperator.DELETION || e.getOperator() == CigarOperator.INSERTION) {
				return true;
			}
		}
		return false;
	}
	/**
	 * Splits a CIGAR list at the largest indel event
	 * @param list cigar
	 * @return
	 */
	public static List<List<CigarElement>> splitAtLargestIndel(List<CigarElement> list) {
		list = decodeNegativeDeletion(list);
		int bestOpCount = 0;
		int bestOffset = list.size();
		int bestLength = 0;
		for (int i = 0; i < list.size(); i++) {
			int size = 0;
			int length = 0;
			for (int j = i; j < list.size(); j++) {
				CigarElement e = list.get(j);
				if (e.getOperator() == CigarOperator.DELETION || e.getOperator() == CigarOperator.INSERTION) {
					size += Math.abs(e.getLength());
					length = j - i + 1;
				} else {
					break;
				}
			}
			if (size > bestOpCount) {
				bestOffset = i;
				bestOpCount = size;
				bestLength = length;
			}
		}
		List<List<CigarElement>> result = new ArrayList<List<CigarElement>>(3);
		result.add(Lists.newArrayList(Iterables.limit(list, bestOffset)));
		result.add(Lists.newArrayList(Iterables.limit(Iterables.skip(list, bestOffset), bestLength)));
		result.add(Lists.newArrayList(Iterables.skip(list, bestOffset + bestLength)));
		if (result.get(2).size() == 0) {
			result.get(2).add(new CigarElement(0, CigarOperator.MATCH_OR_MISMATCH));
		}
		return result;
	}
	public static void addStartSoftClip(List<CigarElement> cigar, int softClippedBaseCount) {
		assert(softClippedBaseCount >= 0);
		while (cigar.get(0).getOperator() == CigarOperator.H) {
			cigar.remove(0);
		}
		if (softClippedBaseCount == 0) return;
		if (cigar.get(0).getOperator() == CigarOperator.SOFT_CLIP) {
			cigar.set(0, new CigarElement(cigar.get(0).getLength() + softClippedBaseCount, CigarOperator.SOFT_CLIP));
		} else {
			cigar.add(0, new CigarElement(softClippedBaseCount, CigarOperator.SOFT_CLIP));
		}
	}
	public static void addEndSoftClip(List<CigarElement> cigar, int softClippedBaseCount) {
		assert(softClippedBaseCount >= 0);
		while (cigar.get(cigar.size() - 1).getOperator() == CigarOperator.H) {
			cigar.remove(cigar.size() - 1);
		}
		if (softClippedBaseCount == 0) return;
		if (cigar.get(cigar.size() - 1).getOperator() == CigarOperator.SOFT_CLIP) {
			cigar.set(cigar.size() - 1, new CigarElement(cigar.get(cigar.size() - 1).getLength() + softClippedBaseCount, CigarOperator.SOFT_CLIP));
		} else {
			cigar.add(new CigarElement(softClippedBaseCount, CigarOperator.SOFT_CLIP));
		}
	}
	public static void trimClipping(List<CigarElement> cigar) {
		while (trimClippingAt(cigar, 0));
		while (trimClippingAt(cigar, cigar.size() - 1));
	}
	private static boolean trimClippingAt(List<CigarElement> cigar, int offset) {
		if (cigar.get(offset).getOperator() == CigarOperator.HARD_CLIP || cigar.get(offset).getOperator() == CigarOperator.SOFT_CLIP) {
			cigar.remove(offset);
			return true;
		}
		return false;
	}
}
