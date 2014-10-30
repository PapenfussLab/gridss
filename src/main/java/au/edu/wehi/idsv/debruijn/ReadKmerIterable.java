package au.edu.wehi.idsv.debruijn;

import htsjdk.samtools.util.SequenceUtil;

import java.util.Iterator;
import java.util.List;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Lists;
import com.google.common.collect.SortedMultiset;
import com.google.common.collect.TreeMultiset;
import com.google.common.primitives.Bytes;

public class ReadKmerIterable implements Iterable<ReadKmer> {
	private final List<Byte> qual;
	private final List<Byte> bases;
	private final int k;
	private final boolean complement;
	public ReadKmerIterable(int k, byte[] bases, byte[] qual) {
		this(k, bases, qual, false, false);
	}
	public ReadKmerIterable(int k, byte[] bases, byte[] qual, boolean reverse, boolean complement) {
		if (bases == null) {
			throw new NullPointerException("Missing read base information");
		}
		if (bases.length < k) {
			throw new IllegalArgumentException("read length must be greater than k");
		}
		this.complement = complement;
		this.k = k;
		this.bases = reverse ? Lists.reverse(Bytes.asList(bases)) : Bytes.asList(bases);
		this.qual = (qual == null || qual.length == 0) ? null : (reverse ? Lists.reverse(Bytes.asList(qual)) : Bytes.asList(qual));
		if (this.qual != null && this.bases.size() != this.qual.size()) {
			throw new IllegalArgumentException("Quality scores information does not match read base information");
		}
	}
	@Override
	public Iterator<ReadKmer> iterator() {
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
		private final SortedMultiset<Byte> basequals = TreeMultiset.<Byte>create();
		@Override
		protected ReadKmer computeNext() {
			if (offset >= bases.size()) return endOfData();
			if (offset == Integer.MIN_VALUE) {
				init();
			} else {
				advance();
			}
			// add 1 to qual to ensure it is always positive
			return new ReadKmer(currentkmer, 1 + (qual == null ? 0 : basequals.firstEntry().getElement()), lastAmbigiousBaseOffset >= offset - k);
		}
		private void advance() {
			currentkmer = KmerEncodingHelper.nextState(k, currentkmer, adjustBase(bases.get(offset)));
			if (qual != null) {
				basequals.remove(qual.get(offset - k));
				basequals.add(qual.get(offset));
			}
			if (KmerEncodingHelper.isAmbiguous(bases.get(offset))) {
				lastAmbigiousBaseOffset = offset;
			}
			offset++;
		}
		private void init() {
			currentkmer = KmerEncodingHelper.picardBaseToEncoded(k, bases);
			if (complement) currentkmer = KmerEncodingHelper.complement(k, currentkmer);
			if (qual != null) {
				for (int i = 0; i < k; i++) {
					basequals.add(qual.get(i));
				}
			}
			offset = k;
			for (int i = 0; i < k; i++) {
				if (KmerEncodingHelper.isAmbiguous(bases.get(i))) {
					lastAmbigiousBaseOffset = i;
				}
			}
		}
	}
}
