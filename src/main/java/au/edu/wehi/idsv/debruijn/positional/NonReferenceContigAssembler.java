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
import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.Models;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.debruijn.DeBruijnGraphBase;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.graph.ScalingHelper;
import au.edu.wehi.idsv.util.IntervalUtil;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;


/**
 * Calls optimal contigs from a positional de Bruijn graph
 * 
 * @author Daniel Cameron
 *
 */
public class NonReferenceContigAssembler extends AbstractIterator<SAMRecordAssemblyEvidence> {
	private Long2ObjectMap<Collection<KmerPathNodeKmerNode>> graphByKmerNode = new Long2ObjectOpenHashMap<Collection<KmerPathNodeKmerNode>>();
	/**
	 * Weights to remove from each kmer path. removeWeight() relies on the following properties of the iterator:
	 * - KmerPathNodes are adjacent, and in ascending order of kmer offset
	 * - traversal does not invoke comparator
	 */ 
	private Object2ObjectRBTreeMap<KmerPathNodeKmerNode, List<KmerNode>> weightToRemove = new Object2ObjectRBTreeMap<KmerPathNodeKmerNode, List<KmerNode>>(KmerPathNodeKmerNode.ByKmerPathNodeOffset);
	private SortedSet<KmerPathNode> graphByPosition = new TreeSet<KmerPathNode>(KmerNodeUtil.ByLastStartEndKmerReferenceWeight);
	private final EvidenceTracker evidenceTracker;
	private final AssemblyEvidenceSource aes;
	private final int maxEvidenceWidth;
	private final int maxAnchorLength;
	private final int k;
	private final int referenceIndex;
	private final KmerPathNodeIteratorInterceptor wrapper;
	public int getReferenceIndex() { return referenceIndex; }
	public NonReferenceContigAssembler(
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
		int targetAnchorLength = Math.max(contig.stream().mapToInt(sn -> sn.length()).sum(), maxAnchorLength);
		KmerPathNodePath startAnchorPath = new KmerPathNodePath(contig.getFirst(), false, targetAnchorLength);
		startAnchorPath.greedyTraverse(true, false);
		ArrayDeque<KmerPathSubnode> startingAnchor = startAnchorPath.headNode().asSubnodes();
		startingAnchor.remove(contig.getFirst());
		while (wrapper.hasNext() && contig.getLast().lastEnd() + targetAnchorLength > wrapper.lastPosition()) {
			// make sure we have enough of the graph loaded so that when
			// we traverse forward, our anchor sequence will be fully defined
			wrapper.next();
		}
		KmerPathNodePath endAnchorPath = new KmerPathNodePath(contig.getLast(), true, targetAnchorLength);
		endAnchorPath.greedyTraverse(true, false);
		ArrayDeque<KmerPathSubnode> endingAnchor = endAnchorPath.headNode().asSubnodes();
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
		List<String> evidenceIds = evidence.stream().map(e -> e.evidenceId()).collect(Collectors.toList());
		SAMRecordAssemblyEvidence assembledContig;
		if (startingAnchor.size() == 0 && endingAnchor.size() == 0) {
			// unanchored
			BreakendSummary be = Models.calculateBreakend(aes.getContext().getLinear(),
					evidence.stream().map(e -> e.breakend()).collect(Collectors.toList()),
					evidence.stream().map(e -> ScalingHelper.toScaledWeight(e.evidenceQuality())).collect(Collectors.toList()));
			assembledContig = AssemblyFactory.createUnanchoredBreakend(aes.getContext(), aes,
					be,
					evidenceIds,
					bases, quals, new int[] { 0, 0 });
		} else if (startingAnchor.size() == 0) {
			// end anchored
			assembledContig = AssemblyFactory.createAnchoredBreakend(aes.getContext(), aes,
					BreakendDirection.Backward, evidenceIds,
					referenceIndex, endingAnchor.getFirst().firstStart(), endingAnchor.stream().mapToInt(n -> n.length()).sum() + k - 1,
					bases, quals, new int[] { 0, 0 });
		} else if (endingAnchor.size() == 0) {
			// start anchored
			assembledContig = AssemblyFactory.createAnchoredBreakend(aes.getContext(), aes,
					BreakendDirection.Forward, evidenceIds,
					referenceIndex, startingAnchor.getLast().lastStart() + k - 1, startingAnchor.stream().mapToInt(n -> n.length()).sum() + k - 1,
					bases, quals, new int[] { 0, 0 });
		} else {
			// left aligned
			assembledContig = AssemblyFactory.createAnchoredBreakpoint(aes.getContext(), aes, evidenceIds,
					referenceIndex, startingAnchor.getLast().lastStart() + k - 1, startingAnchor.stream().mapToInt(n -> n.length()).sum() + k - 1,
					referenceIndex, endingAnchor.getFirst().firstStart(), endingAnchor.stream().mapToInt(n -> n.length()).sum() + k - 1,
					bases, quals, new int[] { 0, 0 });
		}
		// remove all evidence contributing to this assembly from the graph
		removeFromGraph(evidence);
		return assembledContig;
		
	}
	/**
	 * Removes all evidence from the current graph
	 * @param evidence
	 */
	private void removeFromGraph(Set<KmerEvidence> evidence) {
		for (KmerEvidence e : evidence) {
			for (int i = 0; i < e.length(); i++) {
				KmerSupportNode support = e.node(i);
				if (support != null) {
					assert(support.lastEnd() <= wrapper.lastPosition());
					addToRemovalList(support);
				}
			}
		}
		removeWeight();
		if (Defaults.PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS) {
			assert(sanityCheck());
		}
	}
	/**
	 * Remove all indicated weight from the graph
	 */
	private void removeWeight() {
		KmerPathNode current = null;
		List<List<KmerNode>> toRemove = null;
		// WARNING: this loop is dangerous as we're mutating the 
		// comparator value of the current node as we traverse the tree
		// Fortunately, we don't change the tree structure as we traverse
		// and we immediately clear the entire tree once we've traversed
		// so we're ok.
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
		assert(node.length() >= toRemove.size());
		// remove from graph
		removeFromGraph(node);
		for (KmerPathNode split : KmerPathNode.removeWeight(node, toRemove)) {
			addToGraph(split);
		}
	}
	private void addToRemovalList(KmerNode support) {
		// find all KmerPathNodeKmerNode matching our support node
		for (KmerPathNodeKmerNode n : graphByKmerNode.get(support.lastKmer())) {
			if (IntervalUtil.overlapsClosed(support.lastStart(), support.lastEnd(), n.lastStart(), n.lastEnd())) {
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
		Collection<KmerPathNodeKmerNode> list = graphByKmerNode.get(node.lastKmer());
		if (list == null) {
			list = new ArrayList<KmerPathNodeKmerNode>();
			graphByKmerNode.put(node.lastKmer(), list);
		}
		list.add(node);
	}
	private void removeFromGraph(KmerPathNodeKmerNode node) {
		Collection<KmerPathNodeKmerNode> list = graphByKmerNode.get(node.lastKmer());
		if (list == null) return;
		list.remove(node);
		if (list.size() == 0) {
			graphByKmerNode.remove(node.lastKmer());
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
			position = node.firstStart();
			return node;
		}
	}
	public boolean sanityCheck() {
		graphByKmerNode.entrySet().stream().flatMap(e -> e.getValue().stream()).forEach(kn -> { 
			assert(kn.node().isValid());
			assert(graphByPosition.contains(kn.node()));
		});
		for (KmerPathNode n : graphByPosition) {
			assert(n.isValid());
		}
		return true;
	}
}
