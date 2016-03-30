package au.edu.wehi.idsv.debruijn.positional;

import htsjdk.samtools.util.Log;
import it.unimi.dsi.fastutil.longs.Long2ObjectMap;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.longs.LongOpenHashSet;
import it.unimi.dsi.fastutil.longs.LongSet;
import it.unimi.dsi.fastutil.objects.ObjectOpenCustomHashSet;

import java.io.File;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.IdentityHashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
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
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.debruijn.DeBruijnGraphBase;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.graph.ScalingHelper;
import au.edu.wehi.idsv.model.Models;
import au.edu.wehi.idsv.util.IntervalUtil;
import au.edu.wehi.idsv.visualisation.PositionalDeBruijnGraphTracker;
import au.edu.wehi.idsv.visualisation.PositionalExporter;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;


/**
 * Calls optimal contigs from a positional de Bruijn graph
 * 
 * @author Daniel Cameron
 *
 */
public class NonReferenceContigAssembler extends AbstractIterator<SAMRecordAssemblyEvidence> {
	private static final Log log = Log.getInstance(NonReferenceContigAssembler.class);
	/**
	 * Positional size of the loaded subgraph before orphan calling is performed.
	 * This is relatively high as orphaned subgraphs are relatively uncommon
	 * This value is in multiples of maxEvidenceDistance
	 */
	private static final int ORPHAN_EVIDENCE_MULTIPLE = 32;
	private Long2ObjectMap<Collection<KmerPathNodeKmerNode>> graphByKmerNode = new Long2ObjectOpenHashMap<Collection<KmerPathNodeKmerNode>>();
	private SortedSet<KmerPathNode> graphByPosition = new TreeSet<KmerPathNode>(KmerNodeUtil.ByFirstStartKmer);
	private final EvidenceTracker evidenceTracker;
	private final AssemblyEvidenceSource aes;
	private final PathCache caller = new PathCache();
	/**
	 * Worst case scenario is a RP providing single kmer support for contig
	 * read length - (k-1) + max-min fragment size
	 *
	 * ========== contig
	 *          --------- read contributes single kmer to contig  
	 *           \       \  in the earliest position
	 *            \  RP   \
	 *             \       \
	 *              \       \
	 *               ---------
	 *                        ^
	 *                        |
	 * Last position supported by this RP is here. 
	 */
	private final int maxEvidenceDistance;
	private final int maxAnchorLength;
	private final int k;
	private final int referenceIndex;
	private final KmerPathNodeIteratorInterceptor wrapper;
	private long consumed = 0;
	private PositionalDeBruijnGraphTracker exportTracker = null;
	public int getReferenceIndex() { return referenceIndex; }
	/**
	 * Creates a new structural variant positional de Bruijn graph contig assembly for the given chromosome
	 * @param it reads
	 * @param referenceIndex evidence source
	 * @param maxEvidenceDistance maximum distance from the first position of the first kmer of a read,
	 *  to the last position of the last kmer of a read. This should be set to read length plus
	 *  the max-min concordant fragment size
	 * @param maxAnchorLength maximum number of reference-supporting anchor bases to assemble
	 * @param k
	 * @param source assembly source
	 * @param tracker evidence lookup
	 */
	public NonReferenceContigAssembler(
			Iterator<KmerPathNode> it,
			int referenceIndex,
			int maxEvidenceDistance,
			int maxAnchorLength,
			int k,
			AssemblyEvidenceSource source,
			EvidenceTracker tracker) {
		this.wrapper = new KmerPathNodeIteratorInterceptor(it);
		this.maxEvidenceDistance = maxEvidenceDistance;
		this.maxAnchorLength = maxAnchorLength;
		this.k = k;
		this.referenceIndex = referenceIndex;
		this.aes = source;
		this.evidenceTracker = tracker;
	}
	@Override
	protected SAMRecordAssemblyEvidence computeNext() {
		SAMRecordAssemblyEvidence calledContig;
		do {
			removeOrphanedNonReferenceSubgraphs();
			ArrayDeque<KmerPathSubnode> bestContig = findBestContig();
			if (bestContig == null) {
				// no more contigs :(
				if (wrapper.hasNext()) {
					log.error("Sanity check failure: end of contigs called before all evidence loaded");
				}
				if (!graphByPosition.isEmpty()) {
					log.error("Sanity check failure: non-empty graph with no contigs called");
				}
				return endOfData();
			}
			calledContig = callContig(bestContig);
		} while (calledContig == null); // if we filtered out our contig, go back 
		return calledContig;
	}
	private ArrayDeque<KmerPathSubnode> findBestContig() {
		Contig contig = findBestCachedContig();
		if (Defaults.SANITY_CHECK_CONTIG_GENERATION) {
			ContigCaller bestCaller = new BestNonReferenceContigCaller(Iterators.concat(graphByPosition.iterator(), wrapper), maxEvidenceDistance);
			ArrayDeque<KmerPathSubnode> bestContig = bestCaller.bestContig();
			assert(bestContig.toString().equals(contig.toSubnodePath().toString()));
		}
		if (exportTracker != null) {
			exportTracker.trackAssembly(caller, contig);
		}
		return contig.toSubnodePath();
	}
	private int nextPosition() {
		if (!wrapper.hasNext()) return Integer.MAX_VALUE;
		return wrapper.peek().firstStart();
	}
	/**
	 * Finds the best available contig
	 * The underlying graph is loaded until such point that a
	 * contig that is the guaranteed maximal contig supported
	 * by the underlying evidence.
	 * @return best contig containing the given supporting evidence 
	 */
	private Contig findBestCachedContig() {
		Contig bestCachedContig = caller.bestContig(nextPosition(), maxEvidenceDistance);
		while (wrapper.hasNext() && bestCachedContig == null) {
			// by loading all nodes within maxEvidenceDistance, we guarantee
			// that all final nodes of incomplete memoized paths
			// can be fully memoized.
			KmerPathNode n = wrapper.next();
			int loadBefore = n.firstStart() + 2 * maxEvidenceDistance;
			while (wrapper.hasNext() && n.firstStart() <= loadBefore) {
				// The next node gets added to the memoization cache
				// by the iterator intercept wrapper calling addToGraph()
				// so we don't have to do anything except advance the
				// iterator here.
				n = wrapper.next();
			}
			bestCachedContig = caller.bestContig(nextPosition(), maxEvidenceDistance);
		}
		return bestCachedContig;
	}
	private SAMRecordAssemblyEvidence callContig(ArrayDeque<KmerPathSubnode> contig) {
		if (containsKmerRepeat(contig)) {
			// recalculate the called contig, this may break the contig at the repeated kmer
			MisassemblyFixer fixed = new MisassemblyFixer(contig);
			contig = new ArrayDeque<KmerPathSubnode>(fixed.correctMisassignedEvidence(evidenceTracker.support(contig)));
		}
		Set<KmerEvidence> evidence = evidenceTracker.untrack(contig);
		
		int targetAnchorLength = Math.max(contig.stream().mapToInt(sn -> sn.length()).sum(), maxAnchorLength);
		KmerPathNodePath startAnchorPath = new KmerPathNodePath(contig.getFirst(), false, targetAnchorLength + maxEvidenceDistance + contig.getFirst().length());
		startAnchorPath.greedyTraverse(true, false);
		ArrayDeque<KmerPathSubnode> startingAnchor = startAnchorPath.headNode().asSubnodes();
		startingAnchor.removeLast();
		while (wrapper.hasNext() && contig.getLast().lastEnd() + targetAnchorLength + maxEvidenceDistance > wrapper.lastPosition()) {
			// make sure we have enough of the graph loaded so that when
			// we traverse forward, our anchor sequence will be fully defined
			wrapper.next();
		}
		KmerPathNodePath endAnchorPath = new KmerPathNodePath(contig.getLast(), true, targetAnchorLength + maxEvidenceDistance + contig.getLast().length());
		endAnchorPath.greedyTraverse(true, false);
		ArrayDeque<KmerPathSubnode> endingAnchor = endAnchorPath.headNode().asSubnodes();
		endingAnchor.removeFirst();
		
		List<KmerPathSubnode> fullContig = new ArrayList<KmerPathSubnode>(contig.size() + startingAnchor.size() + endingAnchor.size());
		fullContig.addAll(startingAnchor);
		fullContig.addAll(contig);
		fullContig.addAll(endingAnchor);
		
		byte[] bases = KmerEncodingHelper.baseCalls(fullContig.stream().flatMap(sn -> sn.node().pathKmers().stream()).collect(Collectors.toList()), k);
		byte[] quals = DeBruijnGraphBase.kmerWeightsToBaseQuals(k, fullContig.stream().flatMapToInt(sn -> sn.node().pathWeights().stream().mapToInt(Integer::intValue)).toArray());
		assert(quals.length == bases.length);
		// left aligned anchor position although it shouldn't matter since anchoring should be a single base wide
		int startAnchorPosition = startingAnchor.size() == 0 ? 0 : startingAnchor.getLast().lastStart() + k - 1;
		int endAnchorPosition = endingAnchor.size() == 0 ? 0 : endingAnchor.getFirst().firstStart();
		int startAnchorBaseCount = startingAnchor.size() == 0 ? 0 : startingAnchor.stream().mapToInt(n -> n.length()).sum() + k - 1;
		int endAnchorBaseCount = endingAnchor.size() == 0 ? 0 : endingAnchor.stream().mapToInt(n -> n.length()).sum() + k - 1;
		int startBasesToTrim = Math.max(0, startAnchorBaseCount - targetAnchorLength);
		int endingBasesToTrim = Math.max(0, endAnchorBaseCount - targetAnchorLength);
		bases = Arrays.copyOfRange(bases, startBasesToTrim, bases.length - endingBasesToTrim);
		quals = Arrays.copyOfRange(quals, startBasesToTrim, quals.length - endingBasesToTrim);
		
		List<String> evidenceIds = evidence.stream().map(e -> e.evidenceId()).collect(Collectors.toList());
		SAMRecordAssemblyEvidence assembledContig;
		if (startingAnchor.size() == 0 && endingAnchor.size() == 0) {
			assert(startBasesToTrim == 0);
			assert(endingBasesToTrim == 0);
			// unanchored
			BreakendSummary be = Models.calculateBreakend(aes.getContext().getLinear(),
					evidence.stream().map(e -> e.breakend()).collect(Collectors.toList()),
					evidence.stream().map(e -> ScalingHelper.toScaledWeight(e.evidenceQuality())).collect(Collectors.toList()));
			assembledContig = AssemblyFactory.createUnanchoredBreakend(aes.getContext(), aes,
					be,
					evidenceIds,
					bases, quals, new int[] { 0, 0 });
			if (evidence.stream().anyMatch(e -> e.isAnchored())) {
				log.debug(String.format("Unanchored assembly %s contains anchored evidence", assembledContig.getEvidenceID()));
			}
		} else if (startingAnchor.size() == 0) {
			// end anchored
			assembledContig = AssemblyFactory.createAnchoredBreakend(aes.getContext(), aes,
					BreakendDirection.Backward, evidenceIds,
					referenceIndex, endAnchorPosition, endAnchorBaseCount - endingBasesToTrim,
					bases, quals);
		} else if (endingAnchor.size() == 0) {
			// start anchored
			assembledContig = AssemblyFactory.createAnchoredBreakend(aes.getContext(), aes,
					BreakendDirection.Forward, evidenceIds,
					referenceIndex, startAnchorPosition, startAnchorBaseCount - startBasesToTrim,
					bases, quals);
		} else {
			if (startAnchorBaseCount + endAnchorBaseCount >= quals.length) {
				// no unanchored bases - not an SV assembly
				assembledContig = null;
			} else {
				assembledContig = AssemblyFactory.createAnchoredBreakpoint(aes.getContext(), aes, evidenceIds,
						referenceIndex, startAnchorPosition, startAnchorBaseCount - startBasesToTrim,
						referenceIndex, endAnchorPosition, endAnchorBaseCount - endingBasesToTrim,
						bases, quals);
			}
		}
		if (aes.getContext().getConfig().getVisualisation().assembly) {
			try {
				PositionalExporter.exportDot(new File(aes.getContext().getConfig().getVisualisation().directory, "assembly." + assembledContig.getEvidenceID() + ".dot"), k, graphByPosition, fullContig);
			} catch (Exception ex) {
				log.debug(ex, "Error exporting assembly ", assembledContig != null ? assembledContig.getEvidenceID() : "(null)");
			}
		}
		if (aes.getContext().getConfig().getVisualisation().fullSizeAssembly) {
			try {
				PositionalExporter.exportNodeDot(new File(aes.getContext().getConfig().getVisualisation().directory, "assembly.fullsize." + assembledContig.getEvidenceID() + ".dot"), k, graphByPosition, fullContig);
			} catch (Exception ex) {
				log.debug(ex, "Error exporting assembly ", assembledContig != null ? assembledContig.getEvidenceID() : "(null)");
			}
		}
		// remove all evidence contributing to this assembly from the graph
		if (evidence.size() > 0) {
			removeFromGraph(evidence);
		} else {
			log.error(String.format("Sanity check failure: found path with no support. Attempting to recover by direct node removal"));
			for (KmerPathSubnode n : contig) {
				removeFromGraph(n.node());
			}
		}
		return assembledContig;
	}
	private boolean containsKmerRepeat(Collection<KmerPathSubnode> contig) {
		LongSet existing = new LongOpenHashSet();
		for (KmerPathSubnode n : contig) {
			for (int i = 0; i < n.length(); i++) {
				if (!existing.add(n.node().kmer(i))) {
					return true;
				}
			}
			for (long kmer : n.node().collapsedKmers()) {
				if (!existing.add(kmer)) {
					return true;
				}
			}
		}
		return false;
	}
	/**
	 * Removes all evidence from the current graph
	 * @param evidence
	 */
	private void removeFromGraph(Set<KmerEvidence> evidence) {
		assert(!evidence.isEmpty());
		Map<KmerPathNode, List<List<KmerNode>>> toRemove = new IdentityHashMap<KmerPathNode, List<List<KmerNode>>>();
		for (KmerEvidence e : evidence) {
			for (int i = 0; i < e.length(); i++) {
				KmerSupportNode support = e.node(i);
				if (support != null) {
					if (support.lastEnd() > wrapper.lastPosition()) {
						log.error(String.format("Sanity check failure: %s extending to %d removed when input at %d", e, support.lastEnd(), wrapper.lastPosition()));
						// try to recover
					}
					updateRemovalList(toRemove, support);
				}
			}
		}
		Set<KmerPathNode> simplifyCandidates = new ObjectOpenCustomHashSet<KmerPathNode>(new KmerPathNode.HashByFirstKmerStartPositionKmer<KmerPathNode>());
		for (Entry<KmerPathNode, List<List<KmerNode>>> entry : toRemove.entrySet()) {
			removeWeight(entry.getKey(), entry.getValue(), simplifyCandidates);
		}
		simplify(simplifyCandidates);
		if (Defaults.SANITY_CHECK_DE_BRUIJN) {
			assert(sanityCheck());
		}
	}
	/**
	 * Attempts to simplify the given nodes
	 * @param simplifyCandidates
	 */
	private void simplify(Set<KmerPathNode> simplifyCandidates) {
		while (!simplifyCandidates.isEmpty()) {
			simplify(simplifyCandidates.iterator().next(), simplifyCandidates);
		}
	}
	private void simplify(KmerPathNode node, Set<KmerPathNode> simplifyCandidates) {
		simplifyCandidates.remove(node);
		KmerPathNode prev = node.prevToMergeWith();
		if (prev != null && prev.lastEnd() < wrapper.lastPosition()) {
			simplifyCandidates.remove(prev);
			removeFromGraph(node);
			removeFromGraph(prev);
			node.prepend(prev);
			addToGraph(node);
		}
		KmerPathNode next = node.nextToMergeWith();
		if (next != null && next.lastEnd() < wrapper.lastPosition()) {
			simplifyCandidates.remove(next);
			removeFromGraph(node);
			removeFromGraph(next);
			next.prepend(node);
			addToGraph(next);
		}
	}
	private void updateRemovalList(Map<KmerPathNode, List<List<KmerNode>>> toRemove, KmerSupportNode support) {
		Collection<KmerPathNodeKmerNode> kpnknList = graphByKmerNode.get(support.lastKmer());
		if (kpnknList != null) {
			for (KmerPathNodeKmerNode n : kpnknList) {
				if (IntervalUtil.overlapsClosed(support.lastStart(), support.lastEnd(), n.lastStart(), n.lastEnd())) {
					updateRemovalList(toRemove, n, support);
				}
			}
		}
	}
	private void updateRemovalList(Map<KmerPathNode, List<List<KmerNode>>> toRemove, KmerPathNodeKmerNode node, KmerSupportNode support) {
		KmerPathNode pn = node.node();
		List<List<KmerNode>> list = toRemove.get(pn);
		if (list == null) {
			list = new ArrayList<List<KmerNode>>(pn.length());
			toRemove.put(pn, list);
		}
		int offset = node.offsetOfPrimaryKmer();
		while (list.size() <= offset) {
			list.add(null);
		}
		List<KmerNode> evidenceList = list.get(offset); 
		if (evidenceList == null) {
			evidenceList = new ArrayList<KmerNode>();
			list.set(offset, evidenceList);
		}
		evidenceList.add(support);
	}
	private void removeWeight(KmerPathNode node, List<List<KmerNode>> toRemove, Set<KmerPathNode> simplifyCandidates) {
		if (node == null) return;
		assert(node.length() >= toRemove.size());
		simplifyCandidates.addAll(node.next());
		simplifyCandidates.addAll(node.prev());
		simplifyCandidates.remove(node);
		Collection<KmerPathNode> replacementNodes = KmerPathNode.removeWeight(node, toRemove);
		// remove from graph
		removeFromGraph(node);
		for (KmerPathNode split : replacementNodes) {
			if (Defaults.SANITY_CHECK_DE_BRUIJN) {
				assert(evidenceTracker.matchesExpected(new KmerPathSubnode(split)));
			}
			addToGraph(split);
		}
		simplifyCandidates.addAll(replacementNodes);
	}
	private void addToGraph(KmerPathNode node) {
		boolean added = graphByPosition.add(node);
		assert(added);
		for (int i = 0; i < node.length(); i++) {
			addToGraph(new KmerPathNodeKmerNode(node, i));
		}
		for (int i = 0; i < node.collapsedKmers().size(); i++) {
			addToGraph(new KmerPathNodeKmerNode(i, node));
		}
		caller.add(node);
	}
	private void removeFromGraph(KmerPathNode node) {
		caller.remove(node);
		boolean removed = graphByPosition.remove(node);
		assert(removed);
		for (int i = 0; i < node.length(); i++) {
			removeFromGraph(new KmerPathNodeKmerNode(node, i));
		}
		for (int i = 0; i < node.collapsedKmers().size(); i++) {
			removeFromGraph(new KmerPathNodeKmerNode(i, node));
		}
	}
	private void addToGraph(KmerPathNodeKmerNode node) {
		Collection<KmerPathNodeKmerNode> list = graphByKmerNode.get(node.firstKmer());
		if (list == null) {
			list = new ArrayList<KmerPathNodeKmerNode>();
			graphByKmerNode.put(node.firstKmer(), list);
		}
		list.add(node);
	}
	private void removeFromGraph(KmerPathNodeKmerNode node) {
		Collection<KmerPathNodeKmerNode> list = graphByKmerNode.get(node.firstKmer());
		if (list == null) return;
		list.remove(node);
		if (list.size() == 0) {
			graphByKmerNode.remove(node.firstKmer());
		}
	}
	/**
	 * Detects and removes orphaned reference subgraphs.
	 * 
	 * Orphaned reference subgraphs can occur when all non-reference
	 * kmers for a given evidence are merged with reference kmers.
	 * Eg: a sequencing error causes a short soft clip leaf which
	 * is subsequently merged with the reference allele.  
	 * 
	 * No contigs will be called for such evidence since no all
	 * connected kmers are reference kmers.
	 * 
	 * @author Daniel Cameron
	 *
	 */
	private void removeOrphanedNonReferenceSubgraphs() {
		if (graphByPosition.isEmpty()) return;
		if (graphByPosition.first().firstStart() >= wrapper.lastPosition() - ORPHAN_EVIDENCE_MULTIPLE * maxEvidenceDistance) return;
		int lastEnd = Integer.MIN_VALUE;
		List<KmerPathNode> nonRefOrphaned = new ArrayList<KmerPathNode>();
		List<KmerPathNode> nonRefActive = new ArrayList<KmerPathNode>();
		// Since we're calling all non-reference contig
		// all orphaned reference subgraph will eventually
		// have no reference kmers with any overlapping positions
		// This means we don't need to actually calculate the subgraph
		// just wait until we've called all the reference contigs and
		// we're just left the non-reference.
		
		// Note: this approach delays removal until all overlapping non-reference
		// subgraphs do not overlap evidence
		// In actual data this does not occur often so is not too much of an issue.
		for (KmerPathNode n : graphByPosition) {
			if (lastEnd < n.firstStart() - 1) {
				// can't be connected since all active nodes
				// have finished sufficiently early that
				// they can't be connected
				nonRefOrphaned.addAll(nonRefActive);
				nonRefActive.clear();
			}
			if (!n.isReference() || n.lastEnd() >= wrapper.lastPosition() - maxEvidenceDistance) {
				// could connect to a reference node
				nonRefActive.clear();
				break;
			}
			lastEnd = Math.max(lastEnd, n.lastEnd());
			nonRefActive.add(n);
		}
		nonRefOrphaned.addAll(nonRefActive);
		if (!nonRefOrphaned.isEmpty()) {
			Set<KmerEvidence> evidence = evidenceTracker.untrack(nonRefOrphaned.stream().map(n -> new KmerPathSubnode(n)).collect(Collectors.toList()));
			removeFromGraph(evidence);			
			// safety check: did we remove them all?
			for (KmerPathNode n : nonRefOrphaned) {
				if (!n.isValid()) continue;
				if (graphByPosition.contains(n)) {
					log.error("Sanity check failure: %s not removed when clearing orphans (%d evidence found). Attempting to recover by direct node removal", n, evidence.size());
					removeFromGraph(n);
				}
			}
		}
	}
	/**
	 * Intercepts the underlying stream and adds nodes to the graph as they are encountered
	 * @author Daniel Cameron
	 *
	 */
	private class KmerPathNodeIteratorInterceptor implements PeekingIterator<KmerPathNode> {
		private final PeekingIterator<KmerPathNode> underlying;
		private int position;
		public int lastPosition() { return position; }
		public KmerPathNodeIteratorInterceptor(Iterator<KmerPathNode> it) {
			this.underlying = Iterators.peekingIterator(it);
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
			consumed++;
			if (Defaults.SANITY_CHECK_DE_BRUIJN) {
				assert(evidenceTracker.matchesExpected(new KmerPathSubnode(node)));
			}
			addToGraph(node);
			position = node.firstStart();
			return node;
		}
		@Override
		public KmerPathNode peek() {
			return underlying.peek();
		}
		@Override
		public void remove() {
			underlying.remove();
		}
	}
	public boolean sanityCheck() {
		graphByKmerNode.entrySet().stream().flatMap(e -> e.getValue().stream()).forEach(kn -> { 
			assert(kn.node().isValid());
			assert(graphByPosition.contains(kn.node()));
		});
		for (KmerPathNode n : graphByPosition) {
			assert(n.isValid());
			assert(evidenceTracker.matchesExpected(new KmerPathSubnode(n)));
		}
		return true;
	}
	public int tracking_activeNodes() {
		return graphByPosition.size();
	}
	public int tracking_maxKmerActiveNodeCount() {
		return graphByKmerNode.values().stream().mapToInt(x -> x.size()).max().orElse(0);
	}
	public long tracking_underlyingConsumed() {
		return consumed;
	}
	public int tracking_inputPosition() {
		return wrapper.lastPosition();
	}
	public int tracking_firstPosition() {
		if (graphByPosition.size() == 0) return Integer.MAX_VALUE;
		return graphByPosition.first().firstStart();
	}
	public PositionalDeBruijnGraphTracker getExportTracker() {
		return exportTracker;
	}
	public void setExportTracker(PositionalDeBruijnGraphTracker exportTracker) {
		this.exportTracker = exportTracker;
	}
}
