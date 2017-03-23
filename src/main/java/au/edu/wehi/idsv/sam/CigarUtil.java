package au.edu.wehi.idsv.sam;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.lang3.tuple.Pair;

import com.google.common.collect.Lists;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class CigarUtil {
	private CigarUtil() {
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
	 * Number of indels in the given cigar
	 * @param list
	 * @return
	 */
	public static int countIndels(Cigar c) {
		int n = 0;
		for (CigarElement e : c.getCigarElements()) {
			if (e.getOperator() == CigarOperator.DELETION || e.getOperator() == CigarOperator.INSERTION) {
				n++;
			}
		}
		return n;
	}
	/**
	 * Split a cigar immediately after the given read base
	 * @param list cigar
	 * @param readBaseOffset 0-based read offset
	 * @return left and right Cigars for corresponding split read alignment
	 */
	public static Pair<Cigar, Cigar> splitAfterReadPosition(final List<CigarElement> list, final int readBaseOffset) {
		int length = readLength(list);
		if (readBaseOffset >= length) throw new IllegalArgumentException(String.format("Offset %d outside of bounds of read with length %d", readBaseOffset, length));
		List<CigarElement> current = new ArrayList<CigarElement>(list.size());
		List<CigarElement> left = current;
		int elementOffset = readBaseOffset;
		for (CigarElement ce : list) {
			if (ce.getOperator().consumesReadBases() && elementOffset >= 0 && elementOffset < ce.getLength()) {
				// split this element
				current.add(new CigarElement(elementOffset + 1, ce.getOperator()));
				clean(current, false);
				current = new ArrayList<CigarElement>(list.size() - current.size() + 2);
				current.add(new CigarElement(ce.getLength() - elementOffset - 1, ce.getOperator()));
			} else {
				// add to current element
				current.add(ce);
			}
			if (ce.getOperator().consumesReadBases()) {
				elementOffset -= ce.getLength();
			}
		}
		left.add(new CigarElement(length - readBaseOffset - 1, CigarOperator.SOFT_CLIP));
		current.add(0, new CigarElement(readBaseOffset + 1, CigarOperator.SOFT_CLIP));
		clean(left, false);
		clean(current, false);
		return Pair.of(new Cigar(left), new Cigar(current));
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
	public static List<CigarElement> trimClipping(List<CigarElement> cigar) {
		while (trimClippingAt(cigar, 0));
		while (trimClippingAt(cigar, cigar.size() - 1));
		return cigar;
	}
	private static boolean trimClippingAt(List<CigarElement> cigar, int offset) {
		if (cigar.get(offset).getOperator() == CigarOperator.HARD_CLIP || cigar.get(offset).getOperator() == CigarOperator.SOFT_CLIP) {
			cigar.remove(offset);
			return true;
		}
		return false;
	}
	public static int countMappedBases(List<CigarElement> cigar) {
		int mapped = 0;
		for (CigarElement e : cigar) {
			if (e.getOperator().consumesReferenceBases() && e.getOperator().consumesReadBases()) {
				mapped += e.getLength();
			}
		}
		return mapped;
	}
	/**
	 * Counts the number of bases mapped to any reference position in both alignment
	 * @param cigar1
	 * @param cigar2
	 * @return number of bases mapped to reference in both alignments
	 */
	public static int commonReferenceBases(Cigar cigar1, Cigar cigar2) {
		int count = 0;
		int readLength = cigar1.getReadLength();
		assert(readLength == cigar2.getReadLength());
		CigarOperatorIterator it1 = new CigarOperatorIterator(cigar1);
		CigarOperatorIterator it2 = new CigarOperatorIterator(cigar2);
		for (int readBaseOffset = 0; readBaseOffset < readLength; readBaseOffset++) {
			CigarOperator op1;
			CigarOperator op2;
			do {
				op1 = it1.next();
			} while (!op1.consumesReadBases());
			do {
				op2 = it2.next();
			} while (!op2.consumesReadBases());
			if (op1.consumesReferenceBases() && op2.consumesReferenceBases()) {
				count++;
			}
		}
		return count;
	}
	public static int getEndClipLength(List<CigarElement> elements) {
		if (elements == null) return 0;
		int length = 0;
		int i = elements.size() - 1;
		while (i >= 0) {
			switch (elements.get(i).getOperator()) {
			case S:
			case H:
				length += elements.get(i).getLength();	
				i--;
				break;
			default:
				return length;
			}
		}
		return length;
	}
	public static int getEndSoftClipLength(List<CigarElement> elements) {
		if (elements == null) return 0;
		int length = 0;
		int i = elements.size() - 1;
		while (i >= 0) {
			switch (elements.get(i).getOperator()) {
			case S:
				length += elements.get(i).getLength();	
				break;
			case H:
				break;
			default:
				return length;
			}
			i--;
		}
		return length;
	}
	public static int getStartClipLength(List<CigarElement> elements) {
		if (elements == null) return 0;
		int length = 0;
		int i = 0;
		while (i < elements.size()) {
			switch (elements.get(i).getOperator()) {
			case S:
			case H:
				length += elements.get(i).getLength();	
				i++;
				break;
			default:
				return length;
			}
		}
		return length;
	}
	public static int getStartSoftClipLength(List<CigarElement> elements) {
		if (elements == null) return 0;
		int length = 0;
		int i = 0;
		while (i < elements.size()) {
			switch (elements.get(i).getOperator()) {
			case S:
				length += elements.get(i).getLength();
				break;
			case H:
				break;
			default:
				return length;
			}
			i++;
		}
		return length;
	}
	/**
	 * Cleans up a cigar by removing zero width events and merging adjacent elements containing the same operator
	 * @param list cigar to clean
	 * @return copy of cleaned list 
	 */
	public static List<CigarElement> clean(List<CigarElement> list) {
		return clean(list, true);
	}
	/**
	 * Cleans up a cigar by removing zero width events and merging adjacent elements containing the same operator
	 * @param list cigar to clean
	 * @param copy determines whether to copy the list before cleaning
	 * @return copy of cleaned list 
	 */
	public static List<CigarElement> clean(List<CigarElement> list, boolean copy) {
		if (copy) {
			list = new ArrayList<CigarElement>(list);
		}
		removeZeroWidthElements(list);
		convertOpenEndedOperations(list);
		mergeAdjacent(list);
		return list;
	}
	/**
	 * Converts indel operations that do not contain anchoring alignments on both sides to
	 * a soft cliped version of the operator 
	 * @param list
	 */
	private static void convertOpenEndedOperations(List<CigarElement> list) {
		for (int i = list.size() - 1; i >= 0; i--) {
			CigarElement ce = list.get(i);
			if (ce.getOperator() == CigarOperator.SOFT_CLIP || ce.getOperator() == CigarOperator.HARD_CLIP) {
			} else if (ce.getOperator() == CigarOperator.INSERTION) {
				list.set(i, new CigarElement(ce.getLength(), CigarOperator.SOFT_CLIP));
			} else if (ce.getOperator() == CigarOperator.DELETION) {
				list.remove(i);
			} else {
				break;
			}
		}
		for (int i = 0; i < list.size(); i++) {
			CigarElement ce = list.get(i);
			if (ce.getOperator() == CigarOperator.SOFT_CLIP || ce.getOperator() == CigarOperator.HARD_CLIP) {
			} else if (ce.getOperator() == CigarOperator.INSERTION) {
				list.set(i, new CigarElement(ce.getLength(), CigarOperator.SOFT_CLIP));
			} else if (ce.getOperator() == CigarOperator.DELETION) {
				list.remove(i);
				i--;
			} else {
				break;
			}
		}
	}
	private static void removeZeroWidthElements(List<CigarElement> list) {
		for (int i = list.size() - 1; i >= 0; i--) {
			if (list.get(i).getLength() == 0) {
				list.remove(i);
			}
		}
	}
	private static void mergeAdjacent(List<CigarElement> list) {
		for (int i = list.size() - 1; i > 0; i--) {
			if (list.get(i).getOperator() == list.get(i - 1).getOperator()) {
				CigarElement replacement = new CigarElement(list.get(i).getLength() + list.get(i - 1).getLength(), list.get(i).getOperator());
				list.set(i - 1, replacement);
				list.remove(i);
			}
		}
	}
	public static class CigarOperatorIterator implements Iterator<CigarOperator> {
		private Iterator<CigarElement> it;
		private CigarElement currentElement = null;
		private int currentIndex;
		public CigarOperatorIterator(Cigar cigar) {
			this(cigar.getCigarElements());
		}
		public CigarOperatorIterator(List<CigarElement> cigar) {
			it = cigar.iterator();
		}
		@Override
		public boolean hasNext() {
			return it.hasNext() || currentElement != null && currentIndex < currentElement.getLength();
		}
		@Override
		public CigarOperator next() {
			if (currentElement == null || currentIndex >= currentElement.getLength()) {
				currentElement = it.next();
				currentIndex = 0;
			}
			currentIndex++;
			return currentElement.getOperator();
		}
	}
	/**
	 * Removes the given number of read bases from the given cigar
	 * @param cigar CIGAR to trim
	 * @param startCount number of starting reads to trim
	 * @param endCount number of ending reads to trim
	 * @return trimmed cigar
	 */
	public static Cigar trimReadBases(Cigar cigar, int startCount, int endCount) {
		List<CigarElement> cl = Lists.newArrayList(cigar.getCigarElements());
		for (int i = 0; i < cl.size(); i++) {
			CigarElement ce = cl.get(i);
			if (!ce.getOperator().consumesReadBases()) {
				cl.remove(i);
				i--;
			} else {
				if (startCount < ce.getLength()) {
					cl.set(i, new CigarElement(ce.getLength() - startCount, ce.getOperator()));
					break;
				} else {
					cl.remove(i);
					i--;
					startCount -= ce.getLength();
				}
			}
		}
		for (int i = cl.size() - 1; i >= 0; i--) {
			CigarElement ce = cl.get(i);
			if (!ce.getOperator().consumesReadBases()) {
				cl.remove(i);
			} else {
				if (endCount < ce.getLength()) {
					cl.set(i, new CigarElement(ce.getLength() - endCount, ce.getOperator()));
					break;
				} else {
					cl.remove(i);
					endCount -= ce.getLength();
				}
			}
		}
		clean(cl);
		return new Cigar(cl);
	}
	/**
	 * Returns the alignment offset of the given base relative to the starting alignment
	 * @param cigar alignment CIGAR
	 * @param readBaseOffset base offset
	 * @return offset relative to first alignment
	 */
	public static int offsetOf(Cigar cigar, int readBaseOffset) {
		List<CigarElement> cl = cigar.getCigarElements();
		int basesLeft = readBaseOffset;
		int currentAlignmentOffset = 0;
		for (int i = 0; i < cl.size(); i++) {
			CigarElement ce = cl.get(i);
			if (ce.getOperator().consumesReadBases() && ce.getOperator().consumesReferenceBases()) {
				if (basesLeft < ce.getLength()) {
					return currentAlignmentOffset + basesLeft;
				}
			}
			if (ce.getOperator().consumesReferenceBases()) {
				currentAlignmentOffset += ce.getLength();
			}
			if (ce.getOperator().consumesReadBases()) {
				basesLeft = Math.max(0, basesLeft - ce.getLength());
			}
		}
		throw new IllegalArgumentException(String.format("Offset of %d base not defined for %s", readBaseOffset, cigar));
	}
	/**
	 * Converts a gapped alignment into an ungapped sub-alignment 
	 * @param cigar cigar
	 * @param startAnchor consider the first mapped bases as the start of the sub-alignment.
	 * If false, the last mapped based is considered the final base of the sub-alignment
	 * @return cigar containing single contiguous mapped sub-alignment 
	 */
	public static List<CigarElement> asUngapped(Cigar cigar, boolean startAnchor) {
		List<CigarElement> list = new ArrayList<>(cigar.getCigarElements());
		boolean convert = false;
		if (startAnchor) {
			for (int i = 0; i < list.size(); i++) {
				CigarElement ci = list.get(i);
				switch (ci.getOperator()) {
					case H:
					case S:
					case M:
					case X:
					case EQ:
						break;
					case D:
					case I:
					case N:
					case P:
						convert = true;
				}
				if (convert) {
					list.set(i, new CigarElement(ci.getOperator().consumesReadBases() ? ci.getLength() : 0, CigarOperator.SOFT_CLIP));
				}
			}
		} else {
			for (int i = list.size() - 1; i >= 0; i--) {
				CigarElement ci = list.get(i);
				switch (ci.getOperator()) {
					case H:
					case S:
					case M:
					case X:
					case EQ:
						break;
					case D:
					case I:
					case N:
					case P:
						convert = true;
				}
				if (convert) {
					list.set(i, new CigarElement(ci.getOperator().consumesReadBases() ? ci.getLength() : 0, CigarOperator.SOFT_CLIP));
				}
			}
		}
		return clean(list, false);
	}
	/**
	 * Assemblies not including anchoring bases cannot be exactly placed but
	 * still have an expected alignment interval. To handle this, we encode the expected
	 * alignment interval in placeholder anchoring bases. These bases are aligned with the
	 * triplet (1X,*N,1X) such that the read alignment interval matches the expected
	 * breakend position interval. Placeholder Ns are added to the assembly so the aligned
	 * anchoring bases do no consume actual assembled sequence.
	 * 
	 * In the case of intervals of 1 or 2 bases, a single 1X or 2X CIGAR anchor is used.
	 * 
	 * @param cigar
	 * @return width of the interval over which the breakend is expected, or 0 if
	 * the CIGAR does not encode in imprecise breakend assembly.
	 */
	public static int widthOfImprecision(final Cigar cigar) {
		if (cigar == null) return 0;
		// should only have XNX
		int width = 0;
		for (CigarElement ce : cigar.getCigarElements()) {
			switch (ce.getOperator()) {
			case X:
			case N:
				width += ce.getLength();
				break;
			case S:
			case H:
				break;
			default:
				// not in the correct format for an unanchored read
				return 0;
			}
		}
		return width;
	}
	public static int countBases(final Cigar cigar, final CigarOperator operator) {
		if (cigar == null) return 0;
		int n = 0;
		for (CigarElement ce : cigar.getCigarElements()) {
			if (ce.getOperator() == operator) {
				n += ce.getLength();
			}
		}
		return n;
	} 
}
