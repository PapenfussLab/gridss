package au.edu.wehi.idsv.visualisation;

import java.io.File;
import java.util.List;
import java.util.Set;

import au.edu.wehi.idsv.debruijn.DeBruijnPathGraph;
import au.edu.wehi.idsv.graph.PathNode;

public interface DeBruijnPathGraphExporter<T, PN extends PathNode<T>> {

	/**
	 * Takes a snapshot of the path graph at the given time
	 * @param pg
	 */
	public abstract DeBruijnPathGraphExporter<T, PN> snapshot(DeBruijnPathGraph<T, PN> pg);

	public abstract DeBruijnPathGraphExporter<T, PN> contigs(List<List<PN>> assembledContigs);
	
	public abstract DeBruijnPathGraphExporter<T, PN> saveTo(File file);

	public abstract void annotateSubgraphs(List<Set<PN>> subgraphs);

	public abstract void annotateStartingPaths(List<Iterable<PN>> startingPaths);

}