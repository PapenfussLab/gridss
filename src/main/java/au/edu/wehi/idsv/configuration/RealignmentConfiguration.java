package au.edu.wehi.idsv.configuration;

import org.apache.commons.configuration.Configuration;

import au.edu.wehi.idsv.AssemblyEvidence;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.SoftClipEvidence;

import com.google.common.collect.ImmutableList;


public class RealignmentConfiguration {
	public static final String CONFIGURATION_PREFIX = "realignment";
	public RealignmentConfiguration(Configuration config) {
		config = config.subset(CONFIGURATION_PREFIX);
		minLength = config.getInt("minLength");
		minAverageQual = config.getFloat("minAverageQual");
		assemblyIterations = config.getInt("assemblyIterations");
		aligner = config.getString("aligner");
		commandline = ImmutableList.copyOf(config.subset("commandline").getStringArray(aligner));
	}
	/**
	 * Minimum breakend length to be considered for realignment
	 */
	public int minLength;
	/**
	 * Minimum average breakend quality score to be considered for realignment 
	 */
	public float minAverageQual;
	/**
	 * Number of realignment iterations performed
	 */
	public int assemblyIterations = 3;
	public String aligner;
	/**
	 * Aligner command line argument. String substitution assumes parameters are
	 * fastq, reference, thread count
	 * @return Aligner command line to execute
	 */
	public ImmutableList<String> commandline;
	public boolean shouldRealignBreakend(SoftClipEvidence evidence) {
		if (evidence.getBreakendSummary() instanceof BreakpointSummary) return false;
		return evidence.getSoftClipLength() >= minLength
				&& evidence.getAverageClipQuality() >= minAverageQual;
	}
	public boolean shouldRealignBreakend(AssemblyEvidence evidence) {
		if (evidence.getBreakendSummary() instanceof BreakpointSummary) return false;
		return evidence.getBreakendSequence().length >= minLength;
				//&& evidence.getAssemblyQuality() >= minAverageQual;
	}
}
