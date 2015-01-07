package au.edu.wehi.idsv;

import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixLoader;
import jaligner.matrix.MatrixLoaderException;


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
	 * Should assembly progress be tracked
	 */
	public static final boolean VISUALISE_ASSEMBLY_PROGRESS;
	/**
	 * Writes evidence that does not pass filters to intermediate files anyway 
	 */
	//public static final boolean WRITE_FILTERED_EVIDENCE;
	/**
	 * Writes calls that do not pass filters to output files anyway 
	 */
	public static final boolean WRITE_FILTERED_CALLS;
	public static final boolean WRITE_FILTERED_ASSEMBLIES;
	public static final boolean PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS;
	public static final boolean PERFORM_EXPENSIVE_CLIQUE_SANITY_CHECKS;
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
	public static final int READ_PAIR_DOVETAIL_MARGIN;
	/**
	 * Minimum number of mapped bases in any alignment
	 */
	public static final int MIN_BASES_TO_ALIGN;
	/**
	 * Number of bases to extend the reference when realigning assemblies
	 */
	public static final int ASSEMBLY_REALIGNMENT_WINDOW_SIZE;
	public static final Matrix ASSEMBLY_REALIGNMENT_MATRIX;
	public static final float ASSEMBLY_REALIGNMENT_PENALTY_GAPOPEN;
	public static final float ASSEMBLY_REALIGNMENT_PENALTY_GAPEXTEND;
	/**
	 * Maximum subgraph width
	 */
	public static final float MAX_SUBGRAPH_WIDTH_IN_FRAGMENT_SIZE_MULTIPLES;
	public static final int ASYNC_READAHEAD_BUFFER_SIZE;
	public static final int ASYNC_READAHEAD_BUFFERS;
	static {
		VISUALISE_ALL = Boolean.valueOf(System.getProperty("gridss.visualisation.saveall", "false"));
		VISUALISE_TIMEOUTS = Boolean.valueOf(System.getProperty("gridss.visualisation.savetimeouts", "false"));
		VISUALISE_ASSEMBLY_PROGRESS = Boolean.valueOf(System.getProperty("gridss.visualisation.trackAssemblyProgress", "false"));
		//WRITE_FILTERED_EVIDENCE = Boolean.valueOf(System.getProperty("gridss.writeFilteredEvidence", "false"));
		WRITE_FILTERED_CALLS = Boolean.valueOf(System.getProperty("gridss.writeFilteredCalls", "false"));
		WRITE_FILTERED_ASSEMBLIES = Boolean.valueOf(System.getProperty("gridss.writeFilteredAssemblies", "false"));
		PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS = Boolean.valueOf(System.getProperty("gridss.expensiveAsserts.debruijn", "false"));
		PERFORM_EXPENSIVE_CLIQUE_SANITY_CHECKS = Boolean.valueOf(System.getProperty("gridss.expensiveAsserts.clique", "false"));
		COLLAPSE_PATH_MAX_TRAVERSAL = Integer.valueOf(System.getProperty("gridss.debruijn.maxCollapseTraversal", "2097152"));
		BEST_PATH_MAX_TRAVERSAL = Integer.valueOf(System.getProperty("gridss.debruijn.maxPathTraversal", "65536"));
		READ_PAIR_DOVETAIL_MARGIN = Integer.valueOf(System.getProperty("gridss.readpair.dovetailMargin", "4"));
		MIN_BASES_TO_ALIGN = Integer.valueOf(System.getProperty("gridss.readpair.minAlignmentBases", "18"));
		ASSEMBLY_REALIGNMENT_WINDOW_SIZE = Integer.valueOf(System.getProperty("gridss.assembly.realignment.window", "20"));
		String matrixName = System.getProperty("gridss.assembly.realignment.matrix", "EDNAFULL");
		try {
			ASSEMBLY_REALIGNMENT_MATRIX = MatrixLoader.load(matrixName);
		} catch (MatrixLoaderException e) {
			throw new RuntimeException("Unable to load alignment matrix " + matrixName);
		}
		ASSEMBLY_REALIGNMENT_PENALTY_GAPOPEN = Float.valueOf(System.getProperty("gridss.assembly.realignment.gapopen", "14"));
		ASSEMBLY_REALIGNMENT_PENALTY_GAPEXTEND = Float.valueOf(System.getProperty("gridss.assembly.realignment.gapextend", "0.5"));
		MAX_SUBGRAPH_WIDTH_IN_FRAGMENT_SIZE_MULTIPLES = Float.valueOf(System.getProperty("gridss.debruijn.maxSubgraphFragmentWidth", "32"));
		ASYNC_READAHEAD_BUFFERS = Integer.valueOf(System.getProperty("gridss.readahead.buffers", "32"));
		ASYNC_READAHEAD_BUFFER_SIZE = Integer.valueOf(System.getProperty("gridss.readahead.buffersize", "128"));
	}
}