package au.edu.wehi.idsv;

import java.io.File;

import org.apache.commons.lang3.StringUtils;

public class Defaults {
	/**
	 * Assembly visualisation directory
	 */
	public static final File ASSEMBLY_VISUALISATION_DIRECTORY;
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
	static {
		String visDir = System.getProperty("gridss.assemblyVisualisationDirectory");
		if (StringUtils.isNotBlank(visDir)) {
			ASSEMBLY_VISUALISATION_DIRECTORY = new File(visDir);
		} else {
			ASSEMBLY_VISUALISATION_DIRECTORY = null;
		}
		//WRITE_FILTERED_EVIDENCE = Boolean.valueOf(System.getProperty("gridss.writeFilteredEvidence", "false"));
		WRITE_FILTERED_CALLS = Boolean.valueOf(System.getProperty("gridss.writeFilteredCalls", "false"));
		PERFORM_EXPENSIVE_DE_BRUIJN_SANITY_CHECKS = Boolean.valueOf(System.getProperty("gridss.debruijn.expensiveAsserts", "false"));
		COLLAPSE_PATH_MAX_TRAVERSAL = Integer.valueOf(System.getProperty("gridss.debruijn.maxCollapseTraversal", "16777216"));
		BEST_PATH_MAX_TRAVERSAL = Integer.valueOf(System.getProperty("gridss.debruijn.maxPathTraversal", "1048576"));
	}
}