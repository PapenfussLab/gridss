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
	static {
		String visDir = System.getProperty("gridss.assemblyVisualisationDirectory");
		if (StringUtils.isNotBlank(visDir)) {
			ASSEMBLY_VISUALISATION_DIRECTORY = new File(visDir);
		} else {
			ASSEMBLY_VISUALISATION_DIRECTORY = null;
		}
		//WRITE_FILTERED_EVIDENCE = Boolean.valueOf(System.getProperty("gridss.writeFilteredEvidence", "false"));
		WRITE_FILTERED_CALLS = Boolean.valueOf(System.getProperty("gridss.writeFilteredCalls", "false"));
	}
}
