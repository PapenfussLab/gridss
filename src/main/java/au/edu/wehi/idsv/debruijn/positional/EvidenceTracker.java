package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.util.IntervalUtil;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.longs.LongLinkedOpenHashSet;
import it.unimi.dsi.fastutil.longs.LongSortedSet;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;

import java.util.*;
import java.util.stream.Collectors;

import static au.edu.wehi.idsv.Defaults.SANITY_CHECK_EVIDENCE_TRACKER;

/**
 * Tracks evidence provided to a given graph by wrapping a source iterator
 * and tracking evidence emitted by the iterator 
 * 
 * @author Daniel Cameron
 *
 */
public class EvidenceTracker {
	//public static EvidenceTracker TEMP_HACK_CURRENT_TRACKER = null;
	private final Long2ObjectOpenHashMap<LinkedList<KmerSupportNode>> lookup = new Long2ObjectOpenHashMap<>();
	private final Object2ObjectOpenHashMap<String, List<KmerEvidence>> id = new Object2ObjectOpenHashMap<>();
	private long evidenceTotal = 0;
	/**
	 * Tracks evidence emitted from the given iterator
	 */
	public EvidenceTracker() {
	}
	/**
	 * Tracks the given evidence
	 * @param support
	 */
	public KmerSupportNode track(KmerSupportNode support) {
		long kmer = support.lastKmer();
		LinkedList<KmerSupportNode> list = lookup.get(kmer);
		if (list == null) {
			list = new LinkedList<KmerSupportNode>();
			lookup.put(kmer, list);
		}
		list.add(support);
		KmerEvidence ke = support.evidence();
		String evidenceId = ke.evidence().getEvidenceID();
		List<KmerEvidence> idvalue = id.get(evidenceId);
		if (idvalue == null) {
			evidenceTotal++;
			idvalue = new ArrayList<>();
			id.put(evidenceId, idvalue);
		}
		if (!idvalue.contains(ke)) {
			idvalue.add(ke);
		}
		return support;
	}
	/**
	 * Stops tracking all nodes associated with all of the given evidence
	 * @param evidenceSet
	 */
	public Set<KmerEvidence> remove(Set<KmerEvidence> evidenceSet) {
		Set<KmerEvidence> evidenceToRemove = new ObjectOpenHashSet<>();
		LongSortedSet kmersInSet = new LongLinkedOpenHashSet();
		for (KmerEvidence evidence : evidenceSet) {
			addToRemoveList(evidence, evidenceToRemove, kmersInSet);
		}
		for (long kmer : kmersInSet) {
			remove(kmer, evidenceToRemove);
		}
		if (SANITY_CHECK_EVIDENCE_TRACKER) {
			sanityCheck();
		}
		return evidenceToRemove;
	}
	private void addToRemoveList(KmerEvidence evidence, Set<KmerEvidence> removeSet, LongSortedSet kmersInSet) {
		// Need to remove all KmerEvidence associated with the evidence
		// Read pairs can have two: one each of the anchored and unanchored reads
		Collection<KmerEvidence> trackedKmerEvidenceForEvidence = id.remove(evidence.evidence().getEvidenceID());
		if (trackedKmerEvidenceForEvidence == null) {
			// Will happen when we attempt to remove the second KmerEvidence in a read pair
			return;
		}
		removeSet.addAll(trackedKmerEvidenceForEvidence);
		for (KmerEvidence e : trackedKmerEvidenceForEvidence) {
			for (int i = 0; i < e.length(); i++) {
				KmerSupportNode node = e.node(i);
				if (node != null) {
					kmersInSet.add(node.firstKmer());
				}
			}
		}
	}
	/**
	 * Stops tracking all nodes associated with the given evidence 
	 * @param evidence
	 */
	private void remove(long kmer, Collection<KmerEvidence> evidence) {
		LinkedList<KmerSupportNode> list = lookup.get(kmer);
		if (list != null) {
			ListIterator<KmerSupportNode> it = list.listIterator();
			while (it.hasNext()) {
				KmerSupportNode n = it.next();
				if (evidence.contains(n.evidence())) {
					it.remove();
				}
			}
			if (list.size() == 0) {
				lookup.remove(kmer);
			}
		}
	}
	/**
	 * Identifies evidence supporting the given path
	 * @param contig path
	 * @return all evidence supporting the given path
	 */
	public Set<KmerEvidence> support(Collection<KmerPathSubnode> contig) {
		return traverse(contig, false);
	}
	/**
	 * Stops tracking all evidence overlapping the given contig
	 * @param contig contig to stop tracking
	 * @return removed evidence 
	 */
	public Set<KmerEvidence> untrack(Collection<KmerPathSubnode> contig) {
		return traverse(contig, true);
	}
	public Set<KmerEvidence> traverse(Collection<KmerPathSubnode> contig, boolean remove) {
		Set<KmerEvidence> evidence = Defaults.USE_OPTIMISED_ASSEMBLY_DATA_STRUCTURES ? new ObjectOpenHashSet<>() : Collections.newSetFromMap(new IdentityHashMap<KmerEvidence, Boolean>());
		for (KmerPathSubnode sn : contig) {
			int start = sn.firstStart();
			int end = sn.firstEnd();
			for (int i = 0; i < sn.length(); i++) {
				toCollection(evidence, sn.kmer(i), start + i, end + i, remove);
			}
		}
		if (remove) {
			// remove any leftover evidence kmers not on the called path
			// Note: this also picks up the RP mate read when only one read is part of the contig
			evidence = remove(evidence);
		}
		if (SANITY_CHECK_EVIDENCE_TRACKER) {
			sanityCheck();
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
	private void toCollection(Collection<KmerEvidence> collection, long kmer, int start, int end, boolean remove) {
		LinkedList<KmerSupportNode> list = lookup.get(kmer);
		if (list != null) {
			ListIterator<KmerSupportNode> it = list.listIterator();
			while (it.hasNext()) {
				KmerSupportNode n = it.next();
				if (IntervalUtil.overlapsClosed(start, end, n.lastStart(), n.lastEnd())) {
					if (remove) {
						it.remove();
					}
					KmerEvidence e = n.evidence();
					collection.add(e);
				}
			}
		}
	}
	public boolean matchesExpected(KmerPathSubnode pn) {
		for (int i = 0; i < pn.length(); i++) {
			LongArrayList kmers = new LongArrayList();
			kmers.add(pn.kmer(i));
			if (!matchesExpected(pn.weight(i) * pn.width(), kmers, pn.firstStart() + i, pn.firstEnd() + i)) {
				return false;
			}
		}
		return true;
	}
	public boolean matchesExpected(int expectedWidthWeight, LongArrayList kmers, int start, int end) {
		int evidenceWeight = 0;
		for (long kmer : kmers) {
			LinkedList<KmerSupportNode> list = lookup.get(kmer);
			if (list != null) {
				for (KmerSupportNode n : list) {
					evidenceWeight += n.weight() * IntervalUtil.overlapsWidthClosed(start, end, n.lastStart(), n.lastEnd());
				}
			}
		}
		assert(evidenceWeight == expectedWidthWeight);
		return evidenceWeight == expectedWidthWeight;
	}
	public boolean isTracked(String evidenceId) {
		return id.keySet().contains(evidenceId);
	}
	public class PathNodeAssertionInterceptor implements Iterator<KmerPathNode> {
		private final Iterator<KmerPathNode> underlying;
		private final String id;
		public PathNodeAssertionInterceptor(Iterator<KmerPathNode> it, String id) {
			this.underlying = it;
			this.id = id;
		}
		@Override
		public boolean hasNext() {
			return underlying.hasNext();
		}
		@Override
		public KmerPathNode next() {
			KmerPathNode node = underlying.next();
			assert(matchesExpected(new KmerPathSubnode(node)));
			return node;
		}
		@Override
		public String toString() {
			return id;
		}
	}
	public class AggregateNodeAssertionInterceptor implements Iterator<KmerNode> {
		private final Iterator<KmerNode> underlying;
		public AggregateNodeAssertionInterceptor(Iterator<KmerNode> it) {
			this.underlying = it;
		}
		@Override
		public boolean hasNext() {
			return underlying.hasNext();
		}
		@Override
		public KmerNode next() {
			KmerNode node = underlying.next();
			assert(matchesExpected(node.width() * node.weight(), LongArrayList.wrap(new long[] { node.firstKmer() }), node.firstStart(), node.firstEnd()));
			return node;
		}
	}
	public Set<KmerEvidence> getTrackedEvidence() {
		return id.values().stream().flatMap(x -> x.stream()).collect(Collectors.toSet());
	}
	public long tracking_evidenceTotal() {
		return evidenceTotal;
	}
	public long tracking_evidenceActive() {
		return id.size();
	}
	public int tracking_kmerCount() {
		return lookup.size();
	}
	public int tracking_supportNodeCount() {
		return lookup.values().stream().mapToInt(x -> x.size()).sum();
	}
	public int tracking_maxKmerSupportNodesCount() {
		return lookup.values().stream().mapToInt(x -> x.size()).max().orElse(0);
	}
	public void sanityCheck() {
		Set<String> lookupEid = lookup.values()
				.stream()
				.flatMap(ll -> ll.stream())
				.map(ksn -> ksn.evidence().evidence().getEvidenceID())
				.collect(Collectors.toSet());
		Set<String> idEid = id.keySet().stream().collect(Collectors.toSet());
		Set<String> missingInLookup = new HashSet<>(idEid);
		Set<String> missingInIds = new HashSet<>(lookupEid);
		missingInIds.removeAll(idEid);
		missingInLookup.removeAll(lookupEid);
		Set<KmerEvidence> kes = lookup.values()
				.stream()
				.flatMap(ll -> ll.stream())
				.map(ksn -> ksn.evidence())
				.collect(Collectors.toSet());
		List<KmerSupportNode> missingKsn = new ArrayList<>();
		for (KmerEvidence ke : kes) {
			for (int i = 0; i < ke.length(); i++) {
				KmerSupportNode ksn = ke.node(i);
				if (ksn != null) {
					LinkedList<KmerSupportNode> list = lookup.get(ksn.firstKmer());
					if (!list.contains(ksn)) {
						missingKsn.add(ksn);
					}
				}
			}
		}
		if (missingInLookup.size() > 0) {
			throw new IllegalStateException("Missing evidence in lookup");
		}
		if (missingInIds.size() > 0) {
			throw new IllegalStateException("Missing all kmer evidence for evidence in lookup");
		}
		if (missingKsn.size() > 0) {
			throw new IllegalStateException("Partially missing kmer support.");
		}
	}
}
