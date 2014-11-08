package au.edu.wehi.idsv.debruijn;

import htsjdk.samtools.util.SequenceUtil;

import java.util.Collections;
import java.util.Iterator;

import org.apache.commons.lang3.ArrayUtils;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.primitives.UnsignedBytes;

public class ReadKmerIterable implements Iterable<ReadKmer> {
	private final byte[] qual;
	private final byte[] bases;
	private final int k;
	private final boolean complement;
	public ReadKmerIterable(int k, byte[] bases, byte[] qual) {
		this(k, bases, qual, false, false);
	}
	public ReadKmerIterable(int k, byte[] bases, byte[] qual, boolean reverse, boolean complement) {
		if (bases == null) {
			throw new NullPointerException("Missing read base information");
		}
		this.complement = complement;
		this.k = k;
		this.bases = reverse ? reverseCopy(bases) : bases;
		this.qual = (qual == null || qual.length == 0) ? null : (reverse ? reverseCopy(qual) : qual);
		if (this.qual != null && this.bases.length != this.qual.length) {
			throw new IllegalArgumentException("Quality scores information does not match read base information");
		}
	}
	private byte[] reverseCopy(byte[] array) {
		byte[] c = array.clone();
		ArrayUtils.reverse(c);
		return c;
	}
	@Override
	public Iterator<ReadKmer> iterator() {
		if (bases.length < k) return Collections.emptyIterator();
		return new ReadKmerIterator();
	}
	private byte adjustBase(byte base) {
		if (complement) return SequenceUtil.complement(base);
		return base;
	}
	private class ReadKmerIterator extends AbstractIterator<ReadKmer> {
		private int lastAmbigiousBaseOffset = Integer.MIN_VALUE;
		private int offset = Integer.MIN_VALUE;
		private long currentkmer;
		private byte[] baseQualsRotatingBuffer;
		private int rotatingBufferPosition;
		private int minBaseQual;
		@Override
		protected ReadKmer computeNext() {
			if (offset >= bases.length) return endOfData();
			if (offset == Integer.MIN_VALUE) {
				init();
			} else {
				advance();
			}
			// add 1 to qual to ensure it is always positive
			return new ReadKmer(currentkmer, 1 + (qual == null ? 0 : minBaseQual), lastAmbigiousBaseOffset >= offset - k);
		}
		private void advance() {
			currentkmer = KmerEncodingHelper.nextState(k, currentkmer, adjustBase(bases[offset]));
			if (qual != null) {
				addToBuffer(qual[offset]);
			}
			if (KmerEncodingHelper.isAmbiguous(bases[offset])) {
				lastAmbigiousBaseOffset = offset;
			}
			offset++;
		}
		private void addToBuffer(byte qual) {
			assert(minBaseQual <= UnsignedBytes.toInt(baseQualsRotatingBuffer[rotatingBufferPosition]));
			minBaseQual = Math.min(UnsignedBytes.toInt(qual), minBaseQual);
			boolean recalcNeeded = minBaseQual == UnsignedBytes.toInt(baseQualsRotatingBuffer[rotatingBufferPosition]);
			baseQualsRotatingBuffer[rotatingBufferPosition] = qual;
			rotatingBufferPosition = (rotatingBufferPosition + 1) % k;
			if (recalcNeeded) {
				// Only need to recalculate when our current min falls out of the kmer
				recalcBufferMin();
			}
		}
		private void recalcBufferMin() {
			// linear traversal faster than SortedMultiset<Byte> for small k
			minBaseQual = UnsignedBytes.min(baseQualsRotatingBuffer);
		}
		private void init() {
			currentkmer = KmerEncodingHelper.picardBaseToEncoded(k, bases);
			if (complement) currentkmer = KmerEncodingHelper.complement(k, currentkmer);
			if (qual != null) {
				baseQualsRotatingBuffer = new byte[k];
				rotatingBufferPosition = 0;
				System.arraycopy(qual, 0, baseQualsRotatingBuffer, 0, k);
				recalcBufferMin();
			}
			offset = k;
			for (int i = 0; i < k; i++) {
				if (KmerEncodingHelper.isAmbiguous(bases[i])) {
					lastAmbigiousBaseOffset = i;
				}
			}
		}
	}
}
