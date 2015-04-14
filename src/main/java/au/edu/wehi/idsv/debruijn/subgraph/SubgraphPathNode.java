package au.edu.wehi.idsv.debruijn.subgraph;

import java.util.List;

import au.edu.wehi.idsv.debruijn.DeBruijnGraphBase;
import au.edu.wehi.idsv.debruijn.PathNode;

public class SubgraphPathNode extends PathNode<DeBruijnSubgraphNode> {
	private int refCount;
	public SubgraphPathNode(Iterable<Long> path, DeBruijnGraphBase<DeBruijnSubgraphNode> graph) {
		super(path, graph);
	}
	public SubgraphPathNode(SubgraphPathNode unsplit, int startIndex, int length, DeBruijnGraphBase<DeBruijnSubgraphNode> graph) {
		super(unsplit, startIndex, length, graph);
	}
	public SubgraphPathNode(DeBruijnGraphBase<DeBruijnSubgraphNode> graph, Iterable<SubgraphPathNode> path) {
		super(graph, path);
	}
	@Override
	protected void onKmersChanged(DeBruijnGraphBase<DeBruijnSubgraphNode> graph) {
		refCount = 0;
		super.onKmersChanged(graph);
		for (List<Long> kmers : getPathAllKmers()) {
			if (containsReference(graph, kmers)) {
				refCount++;
			}
		}
	}
	public static boolean containsReference(DeBruijnGraphBase<DeBruijnSubgraphNode> graph, List<Long> kmers) {
		for (long kmer : kmers) {
			if (graph.getKmer(kmer).isReference()) return true;
		}
		return false;
	}
	public boolean containsReferenceKmer() { return refCount > 0; }
	public boolean containsNonReferenceKmer() { return refCount < length(); }
	@Override
	protected String printAttributes() {
		return String.format("%s%s", containsReferenceKmer() ? "R" : " ", containsNonReferenceKmer()? "N" : " ") + super.printAttributes();
	}
}
