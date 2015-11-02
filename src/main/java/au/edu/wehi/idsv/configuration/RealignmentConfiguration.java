package au.edu.wehi.idsv.configuration;

import htsjdk.samtools.SAMRecord;

import org.apache.commons.configuration.Configuration;

import au.edu.wehi.idsv.AssemblyEvidence;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.SoftClipEvidence;


public class RealignmentConfiguration {
	public static final String CONFIGURATION_PREFIX = "realignment";
	public RealignmentConfiguration(Configuration config) {
		config = config.subset(CONFIGURATION_PREFIX);
		minLength = config.getInt("minLength");
		minAverageQual = config.getFloat("minAverageQual");
		mapqUniqueThreshold = config.getInt("mapqUniqueThreshold");
		assemblyIterations = config.getInt("assemblyIterations");
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
	 * Minimum MAPQ of realigned segment to be considered uniquely aligned  
	 */
	public int mapqUniqueThreshold;
	/**
	 * Number of realignment iterations performed
	 */
	public int assemblyIterations = 3;
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
	/**
	 * Determines whether the given realignment is considered uniquely mappable
	 * @param realignment realignment record
	 * @return true if the mapping position is unique, false otherwise
	 */
	public boolean realignmentPositionUnique(SAMRecord realignment) {
		return realignment.getMappingQuality() >= mapqUniqueThreshold;
	}
}
