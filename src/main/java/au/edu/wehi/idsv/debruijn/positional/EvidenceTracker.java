package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.longs.LongArrayList;

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Set;

import au.edu.wehi.idsv.util.IntervalUtil;

public class EvidenceTracker {
	private final Long2ObjectOpenHashMap<LinkedList<KmerSupportNode>> lookup = new Long2ObjectOpenHashMap<LinkedList<KmerSupportNode>>();
	/**
	 * Tracks the given evidence
	 * @param evidence
	 */
	public void track(KmerSupportNode support) {
		long kmer = support.kmer();
		LinkedList<KmerSupportNode> list = lookup.get(kmer);
		if (list == null) {
			list = new LinkedList<KmerSupportNode>();
			lookup.put(kmer, list);
		}
		list.add(support);
	}
	/**
	 * Stops tracking all evidence overlapping the given contig
	 * @param contig contig to stop tracking
	 * @return matching evidence
	 */
	public Set<KmerEvidence> untrack(Collection<KmerPathSubnode> contig) {
		HashSet<KmerEvidence> e = new HashSet<KmerEvidence>();
		for (KmerPathSubnode sn : contig) {
			int start = sn.firstKmerStartPosition();
			int end = sn.firstKmerEndPosition();
			for (int i = 0; i < sn.length(); i++) {
				removeToCollection(e, sn.kmer(i), start + i, end + i);
			}
			LongArrayList collapsed = sn.node().collapsedKmers();
			IntArrayList collapsedOffset = sn.node().collapsedKmerOffsets();
			for (int i = 0; i < collapsed.size(); i++) {
				int offset = collapsedOffset.getInt(i);
				removeToCollection(e, collapsed.getLong(i), start + offset, end + offset);
			}
		}
		return e;
	}
	public void removeToCollection(Collection<KmerEvidence> collection, long kmer, int start, int end) {
		LinkedList<KmerSupportNode> list = lookup.get(kmer);
		if (list != null) {
			Iterator<KmerSupportNode> it = list.iterator();
			while (it.hasNext()) {
				KmerSupportNode n = it.next();
				if (IntervalUtil.overlapsClosed(start, end, n.startPosition(), n.endPosition())) {
					it.remove();
					list.add(n);
				}
			}
		}
	}
}
