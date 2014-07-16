package au.edu.wehi.idsv;
/**
 * Counts the number of reads and read pairs providing support for the
 * absence of a structural variation at a given position
 * 
 * @author Daniel Cameron
 *
 */
public interface ReferenceCoverageLookup {

	/**
	 * Number of reference reads providing evidence against a breakend immediately after the given base
	 * @param referenceIndex contig
	 * @param position position immediate before putative breakend
	 * @return number of reads spanning the putative breakend
	 */
	public abstract int readsSupportingNoBreakendAfter(int referenceIndex,
			int position);

	/**
	 * Number of read pairs providing evidence against a breakend immediately after the given base
	 * @param referenceIndex contig
	 * @param position position immediate before putative breakend
	 * @return number of read pairs spanning the putative breakend
	 */
	public abstract int readPairsSupportingNoBreakendAfter(int referenceIndex,
			int position);

}