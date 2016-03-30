package au.edu.wehi.idsv.debruijn.positional;

import java.util.ArrayDeque;
import java.util.Iterator;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

public abstract class ContigCaller {

	/**
	 * Since reference kmers are not scored, calculating 
	 * highest weighted results in a preference for paths
	 * ending at a RP with sequencing errors over a path
	 * anchored to the reference. 
	 * 
	 * To ensure that the anchored paths are scored higher
	 * than the unanchored paths, paths anchored to the
	 * reference are given a score adjustment larger than
	 * the largest expected score.
	 */
	protected static final int ANCHORED_SCORE = Integer.MAX_VALUE >> 2;

	public abstract ArrayDeque<KmerPathSubnode> bestContig();

	protected final PeekingIterator<KmerPathNode> underlying;
	protected final int maxEvidenceWidth;
	public int nextPosition() {
		if (!underlying.hasNext()) return Integer.MAX_VALUE;
		return underlying.peek().firstStart();
	}
	public ContigCaller(
			Iterator<KmerPathNode> it,
			int maxEvidenceWidth) {
		this.underlying = Iterators.peekingIterator(it);
		this.maxEvidenceWidth = maxEvidenceWidth;
	}
}