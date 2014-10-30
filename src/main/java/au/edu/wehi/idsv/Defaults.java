package au.edu.wehi.idsv;


public class Defaults {
	/**
	 * Should assemblies causing timeouts be saved
	 */
	public static final boolean VISUALISE_ALL;
	/**
	 * Should assemblies causing errors be saved
	 */
	public static final boolean VISUALISE_TIMEOUTS;
	/**
	 * Writes evidence that does not pass filters to intermediate files anyway 
	 */
	//public static final boolean WRITE_FILTERED_EVIDENCE;
	/**
	 * Writes calls that do not pass filters to output files anyway 
	 */
	public static final boolean WRITE_FILTERED_CALLS;
	public static final boolean PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS;
	/**
	 * Safety limit to prevent unbounded exponential runtime
	 * when attempt to path collapse highly collected degenerate subgraphs
	 */
	public static final int COLLAPSE_PATH_MAX_TRAVERSAL;
	/**
	 * Safety limit to prevent unbounded exponential runtime
	 * when attempt to path collapse highly collected degenerate subgraphs
	 */
	public static final int BEST_PATH_MAX_TRAVERSAL;
	/**
	 * Maximum subgraph width
	 */
	public static final float MAX_SUBGRAPH_WIDTH_IN_FRAGMENT_SIZE_MULTIPLES;
	static {
		VISUALISE_ALL = Boolean.valueOf(System.getProperty("gridss.visualisation.saveall", "false"));
		VISUALISE_TIMEOUTS = Boolean.valueOf(System.getProperty("gridss.visualisation.savetimeouts", "false"));
		//WRITE_FILTERED_EVIDENCE = Boolean.valueOf(System.getProperty("gridss.writeFilteredEvidence", "false"));
		WRITE_FILTERED_CALLS = Boolean.valueOf(System.getProperty("gridss.writeFilteredCalls", "false"));
		PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS = Boolean.valueOf(System.getProperty("gridss.debruijn.expensiveAsserts", "false"));
		COLLAPSE_PATH_MAX_TRAVERSAL = Integer.valueOf(System.getProperty("gridss.debruijn.maxCollapseTraversal", "16777216"));
		BEST_PATH_MAX_TRAVERSAL = Integer.valueOf(System.getProperty("gridss.debruijn.maxPathTraversal", "1048576"));
		MAX_SUBGRAPH_WIDTH_IN_FRAGMENT_SIZE_MULTIPLES = Float.valueOf(System.getProperty("gridss.debruijn.maxSubgraphFragmentWidth", "32"));
	}
}