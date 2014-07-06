package au.edu.wehi.idsv.debruijn.subgraph;

import htsjdk.samtools.SAMRecord;

import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.PriorityQueue;
import java.util.SortedSet;

import au.edu.wehi.idsv.AssemblyBuilder;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMRecordUtil;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.debruijn.DeBruijnEvidence;
import au.edu.wehi.idsv.debruijn.DeBruijnGraphBase;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.ReadKmer;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

public class DeBruijnReadGraph extends DeBruijnGraphBase<DeBruijnNode> {
	public static final String ASSEMBLER_NAME = "debruijn-s";
	/**
	 * Connected subgraphs
	 */
	private final SortedSet<SubgraphSummary> subgraphs;
	private final int referenceIndex;
	/**
	 * 
	 * @param k
	 * @param direction
	 */
	public DeBruijnReadGraph(ProcessingContext processContext, int referenceIndex, int k, BreakendDirection direction) {
		super(processContext, k, direction);
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
				if (adjKmer == kmer.kmer) continue; // ignore loops back to ourself
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
		for (SubgraphSummary ss : subgraphs) {
			if (ss.getMaxAnchor() < position) {
				SubgraphTraversal subgraph = new SubgraphTraversal(ss.getAnyKmer());
				for (List<Long> contig : subgraph.assembleContigs()) {
					VariantContextDirectedBreakpoint variant = toAssemblyEvidence(contig);
					if (variant != null) {
						contigs.add(variant);
					}
				}
			} else {
				// can break out immediately since subgraphs are sorted by max anchor position
				break;
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
		while (!subgraphs.isEmpty() && subgraphs.first().getMaxAnchor() < position) {
			SubgraphSummary ss = subgraphs.first();
			SubgraphTraversal subgraph = new SubgraphTraversal(ss.getAnyKmer());
			subgraph.remove();
			subgraphs.remove(ss);
		}
	}
	private VariantContextDirectedBreakpoint toAssemblyEvidence(List<Long> contigKmers) {
		int maxsclen = 0;
		int longestSupportingRead = 0;
		int refCount = 0;
		int refAnchor = 0;
		Integer mateAnchor = null;
		// Iterate over reference anchor
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
		// Iterate over breakpoint kmers
		for (long kmer : contigKmers) {
			if (kmers.get(kmer).isReference()) continue;
			Integer mp = kmers.get(kmer).getMatePosition();
			for (SAMRecord support : kmers.get(kmer).getSupportingReads()) {
				longestSupportingRead = Math.max(longestSupportingRead, support.getReadLength());
			}
			if (mateAnchor == null) {
				mateAnchor = mp;
			} else {
				// take closest mate anchor
				mateAnchor = direction == BreakendDirection.Forward ?
						Math.max(mp, mateAnchor) :
						Math.min(mp, mateAnchor);
			}
		}
		AssemblyBuilder builder = debruijnContigAssembly(contigKmers)
				.assemblerName(ASSEMBLER_NAME);
		if (refCount > 0) {
			// anchored read
			builder				
				.referenceAnchor(referenceIndex, refAnchor)
				.maximumSoftClipLength(maxsclen);
		} else if (mateAnchor != null) {
			// inexact breakend
			builder
				.mateAnchor(referenceIndex, mateAnchor)
				.longestSupportingRead(longestSupportingRead);
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
		return builder.makeVariant();
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
			remainingCandidateKmers.remove(start);
			// assemble back until we hit the reference
			for (Long node = greedyPrevState(path.getFirst(), remainingCandidateKmers, null); node != null && !kmers.get(node).isReference(); node = greedyPrevState(node, remainingCandidateKmers, null)) {
				path.addFirst(node);
				remainingCandidateKmers.remove(node); // no longer available for traversal either for us, or for any other contig
			}
			HashSet<Long> referenceVisited = Sets.newHashSet();
			// Extend the contig back along our reference anchor
			for (Long node = greedyPrevState(path.getFirst(), null, referenceVisited); node != null && kmers.get(node).isReference(); node = greedyPrevState(node, null, referenceVisited)) {
				path.addFirst(node);
				referenceVisited.add(node);
			}
			// Then extend forward as far as we can
			for (Long node = greedyNextState(path.getLast(), remainingCandidateKmers, null); node != null && !kmers.get(node).isReference(); node = greedyNextState(node, remainingCandidateKmers, null)) {
				path.addLast(node);
				remainingCandidateKmers.remove(node); // no longer available for traversal either for us, or for any other contig
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
