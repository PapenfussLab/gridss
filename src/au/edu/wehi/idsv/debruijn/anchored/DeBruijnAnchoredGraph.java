package au.edu.wehi.idsv.debruijn.anchored;

import java.util.LinkedList;

import au.edu.wehi.idsv.AssemblyBuilder;
import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.debruijn.DeBruijnNodeBase;
import au.edu.wehi.idsv.debruijn.DeBruijnVariantGraph;
import au.edu.wehi.idsv.debruijn.ReadKmer;
import au.edu.wehi.idsv.debruijn.VariantEvidence;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

public class DeBruijnAnchoredGraph extends DeBruijnVariantGraph<DeBruijnNodeBase> {
	public static final String ASSEMBLER_NAME = "debruijnA";
	private final Multimap<Long, Integer> startkmers = HashMultimap.<Long, Integer>create();
	public DeBruijnAnchoredGraph(ProcessingContext context, AssemblyEvidenceSource source, int k, BreakendDirection direction) {
		super(context, source, k, direction);
	}
	@Override
	protected DeBruijnNodeBase createNode(VariantEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		return new DeBruijnNodeBase(evidence, readKmerOffset, kmer);
	}
	@Override
	protected void onEvidenceAdded(DeBruijnNodeBase graphNode, DeBruijnNodeBase evidenceNode, VariantEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		if (evidence.getReferenceKmerCount() > 0) {
			// we have a fully reference anchored kmer
			if (evidence.isReferenceAnchor(readKmerOffset)) {
				startkmers.put(kmer.kmer, 0);
			}
		} else {
			// no kmer is fully anchored, use the first kmer as that has the most anchored bases
			if (evidence.isFirstKmer(readKmerOffset) && evidence.basesSupportingReference(readKmerOffset) > 0) {
				int offset = getK() - evidence.basesSupportingReference(readKmerOffset);
				startkmers.put(kmer.kmer, offset);
			}
		}
	}
	@Override
	protected void onEvidenceRemoved(DeBruijnNodeBase graphNode, DeBruijnNodeBase evidenceNode, VariantEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		if (evidence.getReferenceKmerCount() > 0) {
			// we have a fully reference anchored kmer
			if (evidence.isReferenceAnchor(readKmerOffset)) {
				startkmers.remove(kmer.kmer, 0);
			}
		} else {
			// no kmer is fully anchored, use the first kmer as that has the most anchored bases
			if (evidence.isFirstKmer(readKmerOffset) && evidence.basesSupportingReference(readKmerOffset) > 0) {
				int offset = getK() - evidence.basesSupportingReference(readKmerOffset);
				startkmers.remove(kmer.kmer, offset);
			}
		}
	}
	public VariantContextDirectedEvidence assembleVariant(int referenceIndex, int position) {
		// debugPrint();
		return greedyAnchoredTraverse(referenceIndex, position);
	}
	/**
	 * Simple greedy traversal starting from the highest weighted starting node
	 * with no repeated nodes
	 * @param position 
	 * @param referenceIndex 
	 * @return
	 */
	private VariantContextDirectedEvidence greedyAnchoredTraverse(int referenceIndex, int position) {
		//debugPrint();
		Long start = bestStartingPosition();
		if (start == null) return null;
		LinkedList<Long> path = greedyTraverse(start);
		VariantContextDirectedEvidence result = pathToAssembly(path, start, referenceIndex, position);
		return result;
	}
	private Long bestStartingPosition() {
		Long best = null;
		long bestScore = -1;
		for (Long state : startkmers.keySet()) {
			long weight = getKmer(state).getWeight();
			if (weight > bestScore) {
				bestScore = weight;
				best = state;
			}
		}
		return best;
	}
	private VariantContextDirectedEvidence pathToAssembly(LinkedList<Long> path, Long breakpointAnchor, int referenceIndex, int position) {
		if (path == null || path.size() == 0) throw new IllegalArgumentException("Invalid path");
		int assemblyLength = path.size() + getK() - 1;
		int offset = getK() - 1;
		int softclipSize = 0;
		for (Long node : path) {
			offset++;
			if (node == breakpointAnchor) {
				softclipSize = assemblyLength - offset;
			}
		}
		if (!startkmers.containsEntry(breakpointAnchor, 0)) {
			// the breakpoint anchor is actually within the kmer, not at the end
			// just grab the the first offset.
			int offsetSize = startkmers.get(breakpointAnchor).iterator().next();
			softclipSize += offsetSize;
		}
		int anchorLen = assemblyLength - softclipSize;
		AssemblyBuilder builder = debruijnContigAssembly(path)
			.assemblerName(ASSEMBLER_NAME)
			.referenceAnchor(referenceIndex, position)
			.anchorLength(anchorLen)
			.contributingEvidence(getSupportingEvidence(path))
			.assembledBaseCount(getSAMRecordBaseCount(path));
		return builder.makeVariant();
	}
	@Override
	public String toString() {
		return super.toString() + String.format("%d start kmers", startkmers.size());
	}
}
