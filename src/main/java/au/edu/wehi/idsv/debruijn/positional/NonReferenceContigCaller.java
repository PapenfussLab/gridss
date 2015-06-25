package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.longs.Long2ObjectMap;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectRBTreeMap;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.AssemblyFactory;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.Models;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.debruijn.DeBruijnGraphBase;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.graph.ScalingHelper;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;


/**
 * Calls optimal contigs from a positional de Bruijn graph
 * 
 * @author Daniel Cameron
 *
 */
public class NonReferenceContigCaller extends AbstractIterator<SAMRecordAssemblyEvidence> {
	private Long2ObjectMap<Collection<KmerPathNodeKmerNode>> graphByKmerNode = new Long2ObjectOpenHashMap<Collection<KmerPathNodeKmerNode>>();
	private Object2ObjectRBTreeMap<KmerPathNodeKmerNode, List<KmerNode>> weightToRemove = new Object2ObjectRBTreeMap<KmerPathNodeKmerNode, List<KmerNode>>(KmerPathNodeKmerNode.ByKmerPathNodeOffset);
	private SortedSet<KmerPathNode> graphByPosition = new TreeSet<KmerPathNode>(KmerPathNode.ByStartEndPositionKmerReferenceWeight);
	private final EvidenceTracker evidenceTracker;
	private final AssemblyEvidenceSource aes;
	private final int maxEvidenceWidth;
	private final int maxAnchorLength;
	private final int k;
	private final int referenceIndex;
	private final KmerPathNodeIteratorInterceptor wrapper;
	public NonReferenceContigCaller(
			Iterator<KmerPathNode> it,
			int referenceIndex,
			int maxEvidenceWidth,
			int maxAnchorLength,
			int k,
			AssemblyEvidenceSource source,
			EvidenceTracker tracker) {
		this.wrapper = new KmerPathNodeIteratorInterceptor(it);
		this.maxEvidenceWidth = maxEvidenceWidth;
		this.maxAnchorLength = maxAnchorLength;
		this.k = k;
		this.referenceIndex = referenceIndex;
		this.aes = source;
		this.evidenceTracker = tracker;
	}
	@Override
	protected SAMRecordAssemblyEvidence computeNext() {
		BestNonReferenceContigCaller bestCaller = new BestNonReferenceContigCaller(Iterators.concat(graphByPosition.iterator(), wrapper), maxEvidenceWidth);
		ArrayDeque<KmerPathSubnode> contig = bestCaller.bestContig();
		if (contig == null) return endOfData();
		return callContig(contig);
	}
	private SAMRecordAssemblyEvidence callContig(ArrayDeque<KmerPathSubnode> contig) {
		ArrayDeque<KmerPathSubnode> startingAnchor = new KmerPathNodePath(contig.getFirst(), false, maxAnchorLength).headNode().asSubnodes();
		startingAnchor.remove(contig.getFirst());
		while (wrapper.hasNext() && contig.getLast().endPosition() + maxAnchorLength > wrapper.lastPosition()) {
			// make sure we have enough of the graph loaded that when
			// we traverse forward, our anchor sequence is fully defined
			wrapper.next();
		}
		ArrayDeque<KmerPathSubnode> endingAnchor = new KmerPathNodePath(contig.getLast(), true, maxAnchorLength).headNode().asSubnodes();
		endingAnchor.remove(contig.getLast());
		
		List<KmerPathSubnode> fullContig = new ArrayList<KmerPathSubnode>(contig.size() + startingAnchor.size() + endingAnchor.size());
		fullContig.addAll(startingAnchor);
		fullContig.addAll(contig);
		fullContig.addAll(endingAnchor);
		
		byte[] bases = KmerEncodingHelper.baseCalls(fullContig.stream().flatMap(sn -> sn.node().pathKmers().stream()).collect(Collectors.toList()), k);
		byte[] quals = DeBruijnGraphBase.kmerWeightsToBaseQuals(k, fullContig.stream().flatMapToInt(sn -> sn.node().pathWeights().stream().mapToInt(Integer::intValue)).toArray());
		assert(quals.length == bases.length);
		// need to rehydrate evidence
		Set<KmerEvidence> evidence = evidenceTracker.untrack(contig);
		removeFromGraph(evidence);
		List<String> evidenceIds = evidence.stream().map(e -> e.evidenceId()).collect(Collectors.toList());
		if (startingAnchor.size() == 0 && endingAnchor.size() == 0) {
			// unanchored
			BreakendSummary be = Models.calculateBreakend(aes.getContext().getLinear(),
					evidence.stream().map(e -> e.breakend()).collect(Collectors.toList()),
					evidence.stream().map(e -> ScalingHelper.toScaledWeight(e.evidenceQuality())).collect(Collectors.toList()));
			return AssemblyFactory.createUnanchoredBreakend(aes.getContext(), aes,
					be,
					evidenceIds,
					bases, quals, new int[] { 0, 0 });
		} else if (startingAnchor.size() == 0) {
			// end anchored
			return AssemblyFactory.createAnchoredBreakend(aes.getContext(), aes,
					BreakendDirection.Backward, evidenceIds,
					referenceIndex, endingAnchor.getFirst().firstKmerStartPosition(), endingAnchor.stream().mapToInt(n -> n.length()).sum(),
					bases, quals, new int[] { 0, 0 });
		} else if (endingAnchor.size() == 0) {
			// start anchored
			return AssemblyFactory.createAnchoredBreakend(aes.getContext(), aes,
					BreakendDirection.Forward, evidenceIds,
					referenceIndex, startingAnchor.getLast().startPosition(), startingAnchor.stream().mapToInt(n -> n.length()).sum(),
					bases, quals, new int[] { 0, 0 });
		} else {
			// left aligned
			return AssemblyFactory.createAnchoredBreakpoint(aes.getContext(), aes, evidenceIds,
					referenceIndex, startingAnchor.getLast().startPosition(), startingAnchor.stream().mapToInt(n -> n.length()).sum(),
					referenceIndex, endingAnchor.getFirst().firstKmerStartPosition(), endingAnchor.stream().mapToInt(n -> n.length()).sum(),
					bases, quals, new int[] { 0, 0 });
		}
		
	}
	/**
	 * Removes all evidence from the current graph
	 * @param evidence
	 */
	private void removeFromGraph(Set<KmerEvidence> evidence) {
		for (KmerEvidence e : evidence) {
			for (int i = 0; i < e.length(); i++) {
				KmerSupportNode support = e.node(i);
				assert(support.endPosition() <= wrapper.lastPosition());
				if (support != null) {
					addToRemovalList(support);
				}
			}
		}
		removeWeight();
	}
	/**
	 * Remove all indicated weight from the graph
	 */
	private void removeWeight() {
		KmerPathNode current = null;
		List<List<KmerNode>> toRemove = null;
		for (Entry<KmerPathNodeKmerNode, List<KmerNode>> entry : weightToRemove.entrySet()) {
			KmerPathNodeKmerNode node = entry.getKey();
			if (node.node() != current) {
				removeWeight(current, toRemove);
				current = node.node();
				toRemove = new ArrayList<List<KmerNode>>(node.node().length());
			}
			int offset = node.offsetOfPrimaryKmer();
			while (toRemove.size() < offset) {
				toRemove.add(null);
			}
			toRemove.add(entry.getValue());
		}
		removeWeight(current, toRemove);
		weightToRemove.clear();
	}
	private void removeWeight(KmerPathNode node, List<List<KmerNode>> toRemove) {
		if (node == null) return;
		// remove from graph
		removeFromGraph(node);
		for (KmerPathNode split : KmerPathNode.removeWeight(node, toRemove)) {
			addToGraph(split);
		}
	}
	private void addToRemovalList(KmerNode support) {
		// find all KmerPathNodeKmerNode matching our support node
		for (KmerPathNodeKmerNode n : graphByKmerNode.get(support.kmer())) {
			if (support.overlaps(n)) {
				List<KmerNode> list = weightToRemove.get(n);
				if (list == null) {
					list = new ArrayList<KmerNode>();
					weightToRemove.put(n, list);
				}
				list.add(support);
			}
		}
	}
	private void addToGraph(KmerPathNode node) {
		graphByPosition.add(node);
		for (int i = 0; i < node.length(); i++) {
			addToGraph(new KmerPathNodeKmerNode(node, i));
		}
		for (int i = 0; i < node.collapsedKmers().size(); i++) {
			addToGraph(new KmerPathNodeKmerNode(i, node));
		}
	}
	private void removeFromGraph(KmerPathNode node) {
		graphByPosition.remove(node);
		for (int i = 0; i < node.length(); i++) {
			removeFromGraph(new KmerPathNodeKmerNode(node, i));
		}
		for (int i = 0; i < node.collapsedKmers().size(); i++) {
			removeFromGraph(new KmerPathNodeKmerNode(i, node));
		}
	}
	private void addToGraph(KmerPathNodeKmerNode node) {
		Collection<KmerPathNodeKmerNode> list = graphByKmerNode.get(node.kmer());
		if (list == null) {
			list = new ArrayList<KmerPathNodeKmerNode>();
			graphByKmerNode.put(node.kmer(), list);
		}
		list.add(node);
	}
	private void removeFromGraph(KmerPathNodeKmerNode node) {
		Collection<KmerPathNodeKmerNode> list = graphByKmerNode.get(node.kmer());
		if (list == null) return;
		list.remove(node);
		if (list.size() == 0) {
			graphByKmerNode.remove(node.kmer());
		}
	}
	/**
	 * Intercepts the underlying stream and adds nodes to the graph as they are encountered
	 * @author cameron.d
	 *
	 */
	private class KmerPathNodeIteratorInterceptor implements Iterator<KmerPathNode> {
		private final Iterator<KmerPathNode> underlying;
		private int position;
		public int lastPosition() { return position; }
		public KmerPathNodeIteratorInterceptor(Iterator<KmerPathNode> it) {
			this.underlying = it;
		}
		@Override
		public boolean hasNext() {
			boolean hasNext = underlying.hasNext();
			if (!hasNext) {
				position = Integer.MAX_VALUE;
			}
			return hasNext;
		}
		@Override
		public KmerPathNode next() {
			KmerPathNode node = underlying.next();
			addToGraph(node);
			position = node.startPosition(0);
			return node;
		}
	}
}
