package au.edu.wehi.idsv.debruijn.subgraph;

import htsjdk.samtools.SAMRecord;

import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;

import au.edu.wehi.idsv.AssemblyBuilder;
import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMRecordUtil;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.debruijn.DeBruijnVariantGraph;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.ReadKmer;
import au.edu.wehi.idsv.debruijn.VariantEvidence;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

public class DeBruijnReadGraph extends DeBruijnVariantGraph<DeBruijnSubgraphNode> {
	public static final String ASSEMBLER_NAME = "debruijn-s";
	/**
	 * Connected subgraphs
	 */
	private final SortedSet<SubgraphSummary> subgraphs;
	private final int referenceIndex;
	private final AssemblyParameters parameters;
	/**
	 * 
	 * @param k
	 * @param direction
	 * @param parameters 
	 */
	public DeBruijnReadGraph(ProcessingContext processContext, AssemblyEvidenceSource source, int referenceIndex, BreakendDirection direction, AssemblyParameters parameters) {
		super(processContext, source, parameters.k, direction);
		this.referenceIndex = referenceIndex;
		this.subgraphs = Sets.newTreeSet(SubgraphSummary.ByMaxAnchor);
		this.parameters = parameters;
	}
	@Override
	protected DeBruijnSubgraphNode createNode(VariantEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		return new DeBruijnSubgraphNode(evidence, readKmerOffset, kmer);
	}
	/**
	 * Sets the subgraph for this new node
	 */
	@Override
	protected void onKmerAdded(long kmer, DeBruijnSubgraphNode node) {
		super.onKmerAdded(kmer, node);
		SubgraphSummary g = null;
		// check if we are connected to an existing subgraph
		for (long adjKmer : KmerEncodingHelper.adjacentStates(getK(), kmer)) {
			if (adjKmer == kmer) continue; // ignore loops back to ourself
			DeBruijnSubgraphNode adjNode = getKmer(adjKmer);
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
			g = new SubgraphSummary(kmer);
		}
		node.setSubgraph(g);
	}
	/**
	 * Updates the containing subgraph details
	 */
	@Override
	protected void onEvidenceAdded(DeBruijnSubgraphNode graphNode, DeBruijnSubgraphNode evidenceNode, VariantEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		SubgraphSummary g = graphNode.getSubgraph();
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
	/**
	 * Assembles contigs which do not have any relevance at or after the given position 
	 * @param position
	 * @return
	 */
	public Iterable<VariantContextDirectedBreakpoint> assembleContigsBefore(int position) {
		List<VariantContextDirectedBreakpoint> contigs = Lists.newArrayList();
		for (SubgraphSummary ss : subgraphs) {
			if (ss.getMaxAnchor() < position) {
				SubgraphPathContigAssembler sca = new SubgraphPathContigAssembler(this, ss.getAnyKmer());
				for (List<Long> contig : sca.assembleContigs(parameters)) {
					VariantContextDirectedBreakpoint variant = toAssemblyEvidence(contig);
					if (variant != null) {
						contigs.add(variant);
					}
				}
			} else {
				// no more subgraphs to process (as subgraphs are sorted by max anchor position)
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
			for (long kmer : reachableFrom(ss.getAnyKmer())) {
				remove(kmer);
			}
			subgraphs.remove(ss);
		}
	}
	/**
	 * Only non-reference reads are considered supporting the breakpoint.
	 */
	@Override
	public Set<SAMRecord> getSupportingSAMRecords(Iterable<Long> path) {
		Set<SAMRecord> reads = Sets.newHashSet();
		for (Long kmer : path) {
			DeBruijnSubgraphNode node = getKmer(kmer);
			if (!node.isReference()) {
				reads.addAll(node.getSupportingReads());
			}
		}
		return reads;
	}
	private VariantContextDirectedBreakpoint toAssemblyEvidence(List<Long> contigKmers) {
		int maxsclen = 0;
		int longestSupportingRead = 0;
		int refCount = 0;
		int refAnchor = 0;
		Integer mateAnchor = null;
		// Iterate over reference anchor
		for (long kmer : contigKmers) {
			if (!getKmer(kmer).isReference()) break;
			refCount++;
			// TODO: choose the best anchor
			refAnchor = getKmer(kmer).getBestReferencePosition();

			// want the max sc len of the final reference kmer 
			maxsclen = 0;
			for (SAMRecord r : getKmer(kmer).getSupportingReads()) {
				maxsclen = Math.max(maxsclen, direction == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipLength(r) : SAMRecordUtil.getStartSoftClipLength(r));
			}
		}
		// Iterate over breakpoint kmers
		for (long kmer : contigKmers) {
			if (getKmer(kmer).isReference()) continue;
			Integer mp = direction == BreakendDirection.Forward ? getKmer(kmer).getMaxMatePosition() : getKmer(kmer).getMinMatePosition();
			for (SAMRecord support : getKmer(kmer).getSupportingReads()) {
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
				.anchorLength(refCount + getK() - 1)
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
}
