package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.longs.LongArrayList;

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.Set;

import au.edu.wehi.idsv.util.IntervalUtil;

/**
 * Tracks evidence provided to a given graph by wrapping a source iterator
 * and tracking evidence emitted by the iterator 
 * 
 * @author cameron.d
 *
 */
public class EvidenceTracker implements Iterator<KmerSupportNode> {
	private final Long2ObjectOpenHashMap<LinkedList<KmerSupportNode>> lookup = new Long2ObjectOpenHashMap<LinkedList<KmerSupportNode>>();
	private final Iterator<KmerSupportNode> underlying;
	/**
	 * Tracks evidence emitted from the given iterator
	 * @param it iterator to track
	 */
	public EvidenceTracker(Iterator<KmerSupportNode> it) {
		this.underlying = it;
	}
	/**
	 * Tracks the given evidence
	 * @param evidence
	 */
	public KmerSupportNode track(KmerSupportNode support) {
		long kmer = support.kmer();
		LinkedList<KmerSupportNode> list = lookup.get(kmer);
		if (list == null) {
			list = new LinkedList<KmerSupportNode>();
			lookup.put(kmer, list);
		}
		list.add(support);
		return support;
	}
	/**
	 * Stops tracking all nodes associated with the given evidence 
	 * @param evidence
	 */
	public void remove(KmerEvidence evidence) {
		for (int i = 0; i < evidence.length(); i++) {
			long kmer = evidence.kmer(i);
			remove(kmer, evidence);
		}
	}
	/**
	 * Stops tracking all nodes associated with the given evidence 
	 * @param evidence
	 */
	public void remove(long kmer, KmerEvidence evidence) {
		LinkedList<KmerSupportNode> list = lookup.get(kmer);
		if (list != null) {
			ListIterator<KmerSupportNode> it = list.listIterator();
			while (it.hasNext()) {
				KmerSupportNode n = it.next();
				if (n.evidence() == evidence) { 
					it.remove();
				}
			}
		}
	}
	/**
	 * Stops tracking all evidence overlapping the given contig
	 * @param contig contig to stop tracking
	 * @return matching evidence
	 */
	public Set<KmerEvidence> untrack(Collection<KmerPathSubnode> contig) {
		HashSet<KmerEvidence> evidence = new HashSet<KmerEvidence>();
		for (KmerPathSubnode sn : contig) {
			int start = sn.firstKmerStartPosition();
			int end = sn.firstKmerEndPosition();
			for (int i = 0; i < sn.length(); i++) {
				removeToCollection(evidence, sn.kmer(i), start + i, end + i);
			}
			LongArrayList collapsed = sn.node().collapsedKmers();
			IntArrayList collapsedOffset = sn.node().collapsedKmerOffsets();
			for (int i = 0; i < collapsed.size(); i++) {
				int offset = collapsedOffset.getInt(i);
				removeToCollection(evidence, collapsed.getLong(i), start + offset, end + offset);
			}
		}
		for (KmerEvidence e : evidence) {
			// remove any leftover evidence kmers not on the called path   
			remove(e);
		}
		return evidence;
	}
	/**
	 * Stops tracking all evidence overlapping the given kmer interval and adds to the given collection
	 * 
	 * @param collection
	 * @param kmer
	 * @param start
	 * @param end
	 */
	public void removeToCollection(Collection<KmerEvidence> collection, long kmer, int start, int end) {
		LinkedList<KmerSupportNode> list = lookup.get(kmer);
		if (list != null) {
			ListIterator<KmerSupportNode> it = list.listIterator();
			while (it.hasNext()) {
				KmerSupportNode n = it.next();
				if (IntervalUtil.overlapsClosed(start, end, n.startPosition(), n.endPosition())) {
					it.remove();
					collection.add(n.evidence());
				}
			}
		}
	}
	@Override
	public boolean hasNext() {
		return underlying.hasNext();
	}
	@Override
	public KmerSupportNode next() {
		// TODO Auto-generated method stub
		return track(underlying.next());
	}
}
