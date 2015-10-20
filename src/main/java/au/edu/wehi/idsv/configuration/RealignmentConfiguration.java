package au.edu.wehi.idsv.configuration;

import au.edu.wehi.idsv.AssemblyEvidence;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.SoftClipEvidence;
import htsjdk.samtools.SAMRecord;


public class RealignmentConfiguration {
	/**
	 * Minimum breakend length to be considered for realignment
	 */
	public int minLength = 25;
	/**
	 * Minimum average breakend quality score to be considered for realignment 
	 */
	public float minAverageQual = new SoftClipConfiguration().minAverageQual;
	/**
	 * Minimum MAPQ of realigned segment to be considered uniquely aligned  
	 */
	public int mapqUniqueThreshold = 10;
	/**
	 * Flag indicating if realignment is required. If true, an error is raised
	 * if a realignment record is required a no record is found.
	 */
	public boolean requireRealignment = true;
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
