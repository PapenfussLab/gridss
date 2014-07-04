package au.edu.wehi.idsv.debruijn.subgraph;

import htsjdk.samtools.SAMRecord;

import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.SortedSet;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.ReadEvidenceAssemblerUtil;
import au.edu.wehi.idsv.SAMRecordUtil;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.VariantContextDirectedBreakpointBuilder;
import au.edu.wehi.idsv.debruijn.DeBruijnEvidence;
import au.edu.wehi.idsv.debruijn.DeBruijnGraphBase;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.ReadKmer;
import au.edu.wehi.idsv.sam.AnomolousReadAssembly;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

public class DeBruijnReadGraph extends DeBruijnGraphBase<DeBruijnNode> {
	public static final String ASSEMBLER_NAME = "debruijnW";
	/**
	 * Connected subgraphs
	 */
	private final SortedSet<SubgraphSummary> subgraphs;
	private final ProcessingContext processContext;
	private final int referenceIndex;
	/**
	 * 
	 * @param k
	 * @param direction
	 */
	public DeBruijnReadGraph(ProcessingContext processContext, int referenceIndex, int k, BreakendDirection direction) {
		super(k, direction);
		this.processContext = processContext;
		this.referenceIndex = referenceIndex;
		this.subgraphs = Sets.newTreeSet(SubgraphSummary.ByMaxAnchor);
	}
	@Override
	protected DeBruijnNode createEmptyNode() {
		return new DeBruijnNode();
	}
	@Override
	protected DeBruijnNode addKmer(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		DeBruijnNode node = super.addKmer(evidence, readKmerOffset, kmer);
		addEvidenceToSubgraph(evidence, addNodeToSubgraph(kmer, node));
		return node;
	}
	private void addEvidenceToSubgraph(DeBruijnEvidence evidence, SubgraphSummary g) {
		int anchorToAdd;
		if (evidence.isDirectlyAnchoredToReference()) {
			anchorToAdd = evidence.getReferenceKmerAnchorPosition();
		} else {
			anchorToAdd = evidence.getMateAnchorPosition();
		}
		if (anchorToAdd < g.getMinAnchor() || anchorToAdd > g.getMaxAnchor()) {
			// we're (possibly) changing the reference bounds of this subgraph
			subgraphs.remove(g);
			g.addAnchor(anchorToAdd);
			subgraphs.add(g);
		} else {
			g.addAnchor(anchorToAdd);
		}
	}
	private SubgraphSummary addNodeToSubgraph(ReadKmer kmer, DeBruijnNode node) {
		SubgraphSummary g = node.getSubgraph();
		if (g == null) {
			// check if we are connected to a subgraph
			for (long adjKmer : KmerEncodingHelper.adjacentStates(k, kmer.kmer)) {
				DeBruijnNode adjNode = kmers.get(adjKmer);
				if (adjNode != null) {
					SubgraphSummary gadj = adjNode.getSubgraph();
					if (g == null) {
						g = gadj;
					} else if (g != gadj) {
						// this kmer merges two subgraphs
						subgraphs.remove(g);
						subgraphs.remove(gadj);
						g = SubgraphSummary.merge(g, gadj);
						subgraphs.add(g);
					}
				}
			}
			if (g == null) {
				// nothing next to us -> new subgraph
				g  = new SubgraphSummary(kmer.kmer);
			}
			node.setSubgraph(g);
		}
		return g;
	}
	/**
	 * Assembles contigs which do not have any relevance at or after the given position 
	 * @param position
	 * @return
	 */
	public Iterable<VariantContextDirectedBreakpoint> assembleContigsBefore(int position) {
		List<VariantContextDirectedBreakpoint> contigs = Lists.newArrayList();
		while (subgraphs.first().getMaxAnchor() < position) {
			SubgraphSummary ss = subgraphs.first();
			SubgraphTraversal subgraph = new SubgraphTraversal(ss.getAnyKmer());
			for (List<Long> contig : subgraph.assembleContigs()) {
				VariantContextDirectedBreakpoint variant = toAssemblyEvidence(contig);
				if (variant) {
					contigs.add();
				}
			}
		}
		Collections.sort(contigs, VariantContextDirectedBreakpoint.ByLocation);
		return contigs;
	}
	/**
	 * Removes all kmers not relevant at or after the given position
	 * @param position
	 */
	public void removeBefore(int position) {
		while (subgraphs.first().getMaxAnchor() < position) {
			SubgraphSummary ss = subgraphs.first();
			SubgraphTraversal subgraph = new SubgraphTraversal(ss.getAnyKmer());
			subgraph.remove();
			subgraphs.remove(ss);
		}
	}
	private VariantContextDirectedBreakpoint toAssemblyEvidence(List<Long> contigKmers) {
		Set<SAMRecord> supportingReads = getSupportingSAMRecords(contigKmers); 
		byte[] bases = getBaseCalls(contigKmers);
		byte[] quals = getBaseQuals(contigKmers);
		int readCount = supportingReads.size();
		int baseCount = getSAMRecordBaseCount(contigKmers);
		
		int maxsclen = 0;
		int longestSupportingRead;
		
		int refCount = 0;
		int refAnchor = 0;
		Integer mateAnchor = null;
		for (long kmer : contigKmers) {
			if (!kmers.get(kmer).isReference()) break;
			refCount++;
			refAnchor = kmers.get(kmer).getReferencePosition();

			// want the max sc len of the final reference kmer 
			maxsclen = 0;
			for (SAMRecord r : kmers.get(kmer).getSupportingReads()) {
				maxsclen = Math.max(maxsclen, direction == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipLength(r) : SAMRecordUtil.getStartSoftClipLength(r));
			}
		}
		for (long kmer : contigKmers) {
			if (kmers.get(kmer).isReference()) continue;
			Integer mp = kmers.get(kmer).getMatePosition();
			if (mateAnchor == null) {
				mp = mateAnchor;
			} else {
				// take closest mate anchor
				mateAnchor = direction == BreakendDirection.Forward ?
						Math.max(mp, mateAnchor) :
						Math.min(mp, mateAnchor);
			}
		}
		AnomolousReadAssembly ara = new AnomolousReadAssembly(ASSEMBLER_NAME, bases, quals, refCount, direction, readCount, baseCount);
		VariantContextDirectedBreakpointBuilder builder;
		if (refCount > 0) {
			builder = ReadEvidenceAssemblerUtil.breakendBuilder(
					processContext,
					ASSEMBLER_NAME,
					referenceIndex,
					refAnchor,
					direction,
					direction == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipBases(ara) : SAMRecordUtil.getStartSoftClipBases(ara),
					direction == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipBaseQualities(ara) : SAMRecordUtil.getStartSoftClipBaseQualities(ara),
					ara.getReadBases(),
					ara.getBaseQualities(),
					ara.getReadCount(),
					ara.getReadBaseCount(),
					ara.getAssemblyBreakpointQuality(),
					maxsclen);
		} else if (mateAnchor != null) {
			// inexact breakend
			builder = ReadEvidenceAssemblerUtil.mateAnchoredBuilder(
					processContext,
					ASSEMBLER_NAME,
					referenceIndex,
					mateAnchor,
					direction,
					ara.getReadBases(),
					ara.getBaseQualities(),
					ara.getReadCount(),
					ara.getReadBaseCount(),
					ara.getAssemblyBreakpointQuality(),
					contigKmers.iterator().next());
		} else {
			// Assembly is neither anchored by breakend nor anchored by mate pair
			// This occurs when the de bruijn graph has a non-reference kmer fork
			// the best path will be taken first which removes the anchor from the
			// forked path. We are left with an unanchored assembly that we can't
			// do anything with
			//                      B-B-B   <- assembly B is unanchored
			//                     /
			//                A-A-A-A-A-A-A-A
			//               /
			// Ref: A-A-A-A-A
			//
			return null;
		}
		VariantContextDirectedBreakpoint variant = builder.make();
		return variant;
	}
	/**
	 * Traverses a given subgraph and iteratively generates contigs from the
	 * putative SVs in a greedy fashion.
	 * 
	 * Each non-reference kmer can be included in only a single assembly.
	 *  
	 * @author Daniel Cameron
	 *
	 */
	private class SubgraphTraversal {
		/**
		 * Remaining non-reference kmers of this subgraph
		 */
		private HashSet<Long> remainingCandidateKmers = Sets.newHashSet();
		/**
		 * Candidate kmers for starting a SV-supporting contig assembly
		 */
		private PriorityQueue<Long> seedCandidates = new PriorityQueue<Long>(256, ByKmerWeight.reverse());
		private final long startKmer;
		private SubgraphTraversal(long startKmer) {
			this.startKmer = startKmer;
		}
		/**
		 * Extracts all SV-supporting contigs from the de bruijn graph
		 * @return kmer contig in subgraph
		 */
		public List<List<Long>> assembleContigs() {
			init();
			List<List<Long>> contigs = Lists.newArrayList();
			for (Long seed = getBestSeed(); seed != null; seed = getBestSeed()) {
				contigs.add(assembleContig(seed));
			}
			return contigs;
		}
		/**
		 * Assembles a contig from the remaining kmers in the subgraph
		 * @param seed starting kmer
		 * @return kmer contig
		 */
		private List<Long> assembleContig(long seed) {
			LinkedList<Long> contigKmers = greedyTraverse(seed);
			return contigKmers;
		}
		/**
		 * Creates a putative SV contig starting from the given kmer.
		 * 
		 * Traversal on non-reference kmers will remove the kmers from the consideration in
		 * further contigs
		 * 
		 * @param start starting kmer seed
		 * @return kmer sequence of contig
		 */
		private LinkedList<Long> greedyTraverse(long start) {
			LinkedList<Long> path = new LinkedList<Long>();
			path.add(start);
			for (Long node = greedyNextState(start, remainingCandidateKmers, null); node != null; node = greedyNextState(node, remainingCandidateKmers, null)) {
				path.addLast(node);
				if (kmers.get(node).isReference()) {
					// We've hit a bubble - time to stop 
					break;
				}
				remainingCandidateKmers.remove(node); // no longer available for traversal either for us, or for any other contig
			}
			HashSet<Long> referenceVisited = Sets.newHashSet();
			boolean inReference = false;
			for (Long node = greedyPrevState(start, remainingCandidateKmers, referenceVisited); node != null; node = greedyPrevState(node, remainingCandidateKmers, referenceVisited)) {
				boolean isReferenceNode = kmers.get(node).isReference(); 
				inReference |= isReferenceNode;
				if (inReference && !isReferenceNode) {
					// Only assemble back to a reference anchor
					// IE: don't traverse back away from the reference once we hit it
					break;
				}
				path.addFirst(node);
				if (inReference) {
					referenceVisited.add(node);
				}
			}
			return path;
		}
		/**
		 * Gets the best kmer to start assembly from
		 * @return highest weighted non-reference kmer, null if no valid starting kmers can be found
		 */
		private Long getBestSeed() {
			while (!seedCandidates.isEmpty()) {
				long kmer = seedCandidates.poll();
				if (remainingCandidateKmers.contains(kmer) && kmers.containsKey(kmer)) {
					return kmer;
				}
			}
			return null;
		}
		/**
		 * Removes the entire subgraph from the de bruijn graph
		 */
		public void remove() {
			HashSet<Long> frontier = Sets.newHashSet();
			frontier.add(startKmer);
			while (!frontier.isEmpty()) {
				long kmer = frontier.iterator().next();
				frontier.remove(kmer);
				// remove kmers from the containing graph
				DeBruijnReadGraph.this.kmers.remove(kmer);
				for (long adjKmer : KmerEncodingHelper.adjacentStates(k, kmer)) {
					if (kmers.containsKey(adjKmer)) {
						frontier.add(adjKmer);
					}
				}
			}
		}
		/**
		 * Initialise the reachability set and best starting kmers for the subgraph
		 */
		private void init() {
			HashSet<Long> frontier = Sets.newHashSet();
			HashSet<Long> referenceVisited = Sets.newHashSet();
			frontier.add(startKmer);
			while (!frontier.isEmpty()) {
				long kmer = frontier.iterator().next();
				frontier.remove(kmer);
				DeBruijnNode node = kmers.get(kmer); 
				if (!node.isReference()) {
					seedCandidates.add(kmer);
					remainingCandidateKmers.add(kmer);
				} else {
					referenceVisited.add(kmer);
				}
				// Add neighbours of this kmer to the frontier
				for (long adjKmer : KmerEncodingHelper.adjacentStates(k, kmer)) {
					// needs to be in the de bruijn graph, and not already processed
					if (kmers.containsKey(adjKmer) && !remainingCandidateKmers.contains(adjKmer) && !referenceVisited.contains(adjKmer)) {
						frontier.add(adjKmer);
					}
				}
			}
		}
	}
}
