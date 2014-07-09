package au.edu.wehi.idsv;

public class AssemblyParameters {
	public AssemblyMethod method = AssemblyMethod.DEBRUIJN_SUBGRAPH;
	/**
	 * De Bruijn graph kmer size
	 */
	public int k = 25;
	/**
	 * Maximum of base mismatches for de bruijn kmer paths to be merged   
	 */
	public int maxBaseMismatchForCollapse = 2;
	/**
	 * Only collapse bubble path with a single entry and exit kmer choice
	 */
	public boolean collapseBubblesOnly = true;
	/**
	 * Allow reuse of reference kmers when assembling subsequent
	 * contigs in an assembly iteration
	 */
	public boolean allReferenceKmerReuse = true;
	/**
	 * Maximum number of contigs per assembly iteration
	 */
	public int maxContigsPerAssembly = 1024;
	/**
	 * Order in which contigs are 
	 * @author Daniel Cameron
	 *
	 */
	public enum ContigAssemblyOrder {
		/**
		 * Greedily assemble traversing to each path containing the highest kmer
		 * 
		 * Starting kmer is the maximum weight non-reference kmer 
		 */
		GreedyMaxKmer,
		/**
		 * Single path with maximal kmer weight
		 */
		OptimalMaxKmer,
	}
	public ContigAssemblyOrder assemblyOrder = ContigAssemblyOrder.GreedyMaxKmer;
}
