package au.edu.wehi.idsv.util;

import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap.Entry;
import it.unimi.dsi.fastutil.ints.Int2DoubleRBTreeMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleSortedMap;
import it.unimi.dsi.fastutil.ints.Int2IntLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.IntBidirectionalIterator;
import it.unimi.dsi.fastutil.ints.IntSortedSet;
import it.unimi.dsi.fastutil.objects.ObjectBidirectionalIterator;

/**
 * Accumulates weighted values according to a predefined set of interval bins
 * @author Daniel Cameron
 *
 */
public class IntervalAccumulator {
	private final Int2DoubleSortedMap bins;
	private Int2IntMap binWidth;
	private final int firstBinStart;
	private final int lastBinEnd;
	public IntervalAccumulator(int start, int end, int binSize) {
		int[] binStart = new int[(end - start) / binSize + 1];
		for (int i = 0; i < binStart.length; i++) {
			binStart[i] = binSize * i + start;
		}
		double[] value = new double[binStart.length];
		this.bins = new Int2DoubleRBTreeMap(binStart, value);
		this.firstBinStart = start;
		this.lastBinEnd = end;
	}
	/**
	 * Splits bins at the given 1-based position such that
	 * a new bin starts at the given position.
	 * 
	 * Warning: must be called before add() as it resets the 
	 * count for the new bin to zero.
	 */
	public void splitBin(int position) {
		if (binWidth != null) {
			throw new IllegalStateException("Must be called before finaliseBins()");
		}
		bins.put(position, 0);
	}
	/**
	 * Indicates that the bins have been finalised and no more changes will occur
	 */
	public void finaliseBins() {
		// precompute bin widths
		binWidth = new Int2IntLinkedOpenHashMap(bins.size());
		IntBidirectionalIterator it = bins.keySet().iterator();
		int lastStart = it.nextInt();
		while (it.hasNext()) {
			int start = it.nextInt();
			binWidth.put(lastStart, start - lastStart);
			lastStart = start;
		}
		binWidth.put(lastStart, lastBinEnd - lastStart + 1);
	}
	/**
	 * Add the given value at all positions in the given interval
	 * @param start start position (inclusive)
	 * @param end end position (inclusive)
	 * @param value value to add
	 */
	public void add(int start, int end, double value) {
		if (binWidth == null) {
			throw new IllegalStateException("Must be called after finaliseBins()");
		}
		if (end < start) {
			throw new IllegalArgumentException("end cannot be before start");
		}
		//double perBaseValue = value / (end - start + 1);
		start = Math.max(start, firstBinStart);
		end = Math.min(end, lastBinEnd);
		int startBinStartKey = bins.headMap(start + 1).lastIntKey();
		int endBinStartKey = bins.headMap(end + 1).lastIntKey();
		Int2DoubleSortedMap toUpdate = bins.subMap(startBinStartKey, endBinStartKey + 1);
		ObjectBidirectionalIterator<Int2DoubleMap.Entry> it = toUpdate.int2DoubleEntrySet().iterator();
		int widthLeft = end - start + 1;
		while (widthLeft > 0) {
			assert(it.hasNext());
			Int2DoubleMap.Entry entry = it.next();
			int entryStartPosition = entry.getIntKey();
			int overlap = IntervalUtil.overlapsWidthClosed(start, end, entryStartPosition, entryStartPosition + getBinSize(entryStartPosition) - 1);
			entry.setValue(entry.getDoubleValue() + overlap * value);
			widthLeft -= overlap;
		}
	}
	public int getBinSize(int binStart) {
		return binWidth.get(binStart);
	}
	/**
	 * Gets the average value for the given bin
	 * @param binStart
	 * @return
	 */
	public double getMeanValue(int binStart) {
		return bins.get(binStart) / getBinSize(binStart);
	}
	public IntSortedSet getBinStarts() {
		return bins.keySet();
	}
	public ObjectBidirectionalIterator<Entry> iterator() {
		return bins.int2DoubleEntrySet().iterator();
	}
}