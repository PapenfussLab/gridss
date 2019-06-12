package au.edu.wehi.idsv.visualisation;

import au.edu.wehi.idsv.debruijn.DeBruijnNodeBase;

import java.io.File;

/**
 * Exports the generated de Bruijn graphs for visualisation
 * @author Daniel Cameron
 *
 * @param <T>
 */
public interface DeBruijnGraphExporter<T extends DeBruijnNodeBase> {
	/*
	 * Processing step time
	 */
	public abstract DeBruijnGraphExporter<T> setTime(int currentTime);

	public abstract DeBruijnGraphExporter<T> trackChanges(long kmer, T kmerNode);

	public abstract DeBruijnGraphExporter<T> saveTo(File file);

}