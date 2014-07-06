package au.edu.wehi.idsv.debruijn.anchored;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;

import au.edu.wehi.idsv.AssemblyBuilder;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.debruijn.DeBruijnEvidence;
import au.edu.wehi.idsv.debruijn.DeBruijnGraphBase;
import au.edu.wehi.idsv.debruijn.DeBruijnNodeBase;
import au.edu.wehi.idsv.debruijn.ReadKmer;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

public class DeBruijnReadGraph extends DeBruijnGraphBase<DeBruijnNodeBase> {
	public static final String ASSEMBLER_NAME = "debruijnA";
	private final Multimap<Long, Integer> startkmers = HashMultimap.<Long, Integer>create();
	public DeBruijnReadGraph(ProcessingContext context, int k, BreakendDirection direction) {
		super(context, k, direction);
	}
	@Override
	protected DeBruijnNodeBase createEmptyNode() {
		return new DeBruijnNodeBase();
	}
	@Override
	protected DeBruijnNodeBase addKmer(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		DeBruijnNodeBase node = super.addKmer(evidence, readKmerOffset, kmer);
		if (evidence.getReferenceKmerCount() > 0) {
			// we have a fully reference anchored kmer
			if (evidence.isReferenceAnchor(readKmerOffset)) {
				startkmers.put(kmer.kmer, 0);
			}
		} else {
			// no kmer is fully anchored, use the first kmer as that has the most anchored bases
			if (evidence.isFirstKmer(readKmerOffset) && evidence.basesSupportingReference(readKmerOffset) > 0) {
				int offset = k - evidence.basesSupportingReference(readKmerOffset);
				startkmers.put(kmer.kmer, offset);
			}
		}
		return node;
	}
	@Override
	protected DeBruijnNodeBase removeKmer(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		DeBruijnNodeBase node = super.removeKmer(evidence, readKmerOffset, kmer);
		if (evidence.getReferenceKmerCount() > 0) {
			// we have a fully reference anchored kmer
			if (evidence.isReferenceAnchor(readKmerOffset)) {
				startkmers.remove(kmer.kmer, 0);
			}
		} else {
			// no kmer is fully anchored, use the first kmer as that has the most anchored bases
			if (evidence.isFirstKmer(readKmerOffset) && evidence.basesSupportingReference(readKmerOffset) > 0) {
				int offset = k - evidence.basesSupportingReference(readKmerOffset);
				startkmers.remove(kmer.kmer, offset);
			}
		}
		return node;
	}
	public VariantContextDirectedBreakpoint assembleVariant(int referenceIndex, int position) {
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
	private VariantContextDirectedBreakpoint greedyAnchoredTraverse(int referenceIndex, int position) {
		//debugPrint();
		Long start = bestStartingPosition();
		if (start == null) return null;
		LinkedList<Long> path = greedyTraverse(start);
		VariantContextDirectedBreakpoint result = pathToAssembly(path, start, referenceIndex, position);
		return result;
	}
	private Long bestStartingPosition() {
		Long best = null;
		long bestScore = -1;
		for (Long state : startkmers.keySet()) {
			long weight = kmers.get(state).getWeight();
			if (weight > bestScore) {
				bestScore = weight;
				best = state;
			}
		}
		return best;
	}
	private VariantContextDirectedBreakpoint pathToAssembly(LinkedList<Long> path, Long breakpointAnchor, int referenceIndex, int position) {
		if (path == null || path.size() == 0) throw new IllegalArgumentException("Invalid path");
		int assemblyLength = path.size() + k - 1;
		int offset = k - 1;
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
			.maximumSoftClipLength(getMaxSoftClipLength(path, anchorLen));
		return builder.makeVariant();
	}
	private LinkedList<Long> greedyTraverse(Long start) {
		LinkedList<Long> path = new LinkedList<Long>();
		Set<Long> visited = new HashSet<Long>();
		path.add(start);
		visited.add(start);
		for (Long node = greedyPrevState(start, null, visited); node != null; node = greedyPrevState(node, null, visited)) {
			path.addFirst(node);
			visited.add(node);
		}
		for (Long node = greedyNextState(start, null, visited); node != null; node = greedyNextState(node, null, visited)) {
			path.addLast(node);
			visited.add(node);
		}
		return path;
	}
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder(super.toString());
		sb.append(String.format("%d start kmers", startkmers.size()));
		return sb.toString();
	}
}
