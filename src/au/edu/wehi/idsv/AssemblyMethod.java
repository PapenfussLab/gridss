package au.edu.wehi.idsv;

public enum AssemblyMethod {
	/**
	 * Independently assemblies a De Bruijn for every genomic position
	 */
	DEBRUIJN_PER_POSITION,
	/**
	 * Assembles local De Bruijn subgraphs, only allowing each kmer to
	 * be used once.  
	 */
	DEBRUIJN_SUBGRAPH,
}
