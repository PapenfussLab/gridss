package au.edu.wehi.idsv.debruijn.positional;

import java.io.File;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.Collection;
import java.util.Iterator;

import com.google.common.collect.Iterators;


/**
 * Wrapper for BestNonReferenceContigCaller allow a single
 * instance to call multiple contigs
 * 
 * @author Daniel Cameron
 *
 */
public class SingleContigCaller extends ContigCaller {
	private final Collection<KmerPathNode> loaded;
	private BestNonReferenceContigCaller caller = null;
	public SingleContigCaller(
			Collection<KmerPathNode> loaded,
			Iterator<KmerPathNode> it,
			int maxEvidenceWidth) {
		super(it, maxEvidenceWidth);
		this.loaded = loaded;
	}
	@Override
	public ArrayDeque<KmerPathSubnode> bestContig() {
		caller = new BestNonReferenceContigCaller(Iterators.concat(loaded.iterator(), underlying), maxEvidenceWidth);
		return caller.bestContig();
	}
	public int tracking_contigCount() {
		return caller.tracking_contigCount();
	}
	public int tracking_contigFirstPosition() {
		return caller.tracking_contigFirstPosition();
	}
	public long tracking_underlyingConsumed() {
		return caller.tracking_underlyingConsumed();
	}
	public int tracking_memoizedNodeCount() {
		return caller.tracking_memoizedNodeCount();
	}
	public int tracking_frontierSize() {
		return caller.tracking_frontierSize();
	}
	public int tracking_unprocessedStartNodeCount() {
		return caller.tracking_unprocessedStartNodeCount();
	}
	@Override
	public void exportState(File file) throws IOException {
		if (caller != null) {
			caller.exportState(file);
		}
	}
}
