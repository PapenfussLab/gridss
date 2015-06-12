package au.edu.wehi.idsv.debruijn.positional;

import java.util.Iterator;
import java.util.PriorityQueue;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

public class SupportNodeIterator implements Iterator<KmerSupportNode> {
	private final PeekingIterator<Evidence> underlying;
	private final int maxSupportWidth;
	private final PriorityQueue<KmerSupportNode> buffer = new PriorityQueue<KmerSupportNode>(1024, KmerNode.ByStartPosition);
	public SupportNodeIterator(Iterator<Evidence> it, int maxSupportWidth) {
		this.underlying = Iterators.peekingIterator(it);
		this.maxSupportWidth = maxSupportWidth;
	}
	private void process(Evidence e) {
		for (int i = 0; i < e.length(); i++) {
			buffer.add(e.node(i));
		}
	}
	@Override
	public boolean hasNext() {
		ensureBuffer();
		return !buffer.isEmpty();
	}
	@Override
	public KmerSupportNode next() {
		ensureBuffer();
		return buffer.poll();
	}
	private void ensureBuffer() {
		while (underlying.hasNext() && buffer.peek().startPosition() + maxSupportWidth >= underlying.peek().startPosition()) {
			process(underlying.next());
		}
	}
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
}