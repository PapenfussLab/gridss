package au.edu.wehi.socrates.debruijn;

import java.util.Iterator;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.SortedMultiset;
import com.google.common.collect.TreeMultiset;

public class ReadKmerIterable implements Iterable<ReadKmer> {
	private final byte[] qual;
	private final byte[] bases;
	private final int k;
	public ReadKmerIterable(int k, byte[] bases, byte[] qual) {
		this.k = k;
		this.bases = bases;
		this.qual = (qual == null || qual.length == 0) ? null : qual;
		if (this.bases == null) {
			throw new NullPointerException("Missing read base information");
		}
		if (this.bases.length < k) {
			throw new IllegalArgumentException("read length must be greater than k");
		}
		if (this.qual != null && this.bases.length != this.qual.length) {
			throw new IllegalArgumentException("Quality scores information does not match read base information");
		}
	}
	@Override
	public Iterator<ReadKmer> iterator() {
		return new ForwardIterator();
	}
	private class ForwardIterator extends AbstractIterator<ReadKmer> {
		private int offset = -1;
		private long currentkmer;
		private final SortedMultiset<Byte> basequals = TreeMultiset.<Byte>create();
		@Override
		protected ReadKmer computeNext() {
			if (offset >= bases.length) return endOfData();
			if (offset < 0) {
				init();
			} else {
				advance();
			}
			// add 1 to qual to ensure it is always positive
			return new ReadKmer(currentkmer, 1 + (qual == null ? 0 : basequals.firstEntry().getElement()));
		}
		private void advance() {
			currentkmer = KmerEncodingHelper.nextState(k, currentkmer, bases[offset]);
			if (qual != null) {
				basequals.remove(qual[offset - k]);
				basequals.add(qual[offset]);
			}
			offset++;
		}
		private void init() {
			currentkmer = KmerEncodingHelper.picardBaseToEncoded(k, bases);
			if (qual != null) {
				for (int i = 0; i < k; i++) {
					basequals.add(qual[i]);
				}
			}
			offset = k;
		}
	}
}
