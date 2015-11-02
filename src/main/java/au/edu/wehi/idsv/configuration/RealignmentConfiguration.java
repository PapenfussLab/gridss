package au.edu.wehi.idsv.configuration;

import java.util.List;

import org.apache.commons.configuration.Configuration;

import com.google.common.collect.ImmutableList;

import au.edu.wehi.idsv.AssemblyEvidence;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.SoftClipEvidence;
import htsjdk.samtools.SAMRecord;


public class RealignmentConfiguration {
	public static final String CONFIGURATION_PREFIX = "realignment";
	public RealignmentConfiguration(Configuration config) {
		config = config.subset(CONFIGURATION_PREFIX);
		minLength = config.getInt("minLength");
		minAverageQual = config.getFloat("minAverageQual");
		mapqUniqueThreshold = config.getInt("mapqUniqueThreshold");
		assemblyIterations = config.getInt("assemblyIterations");
		aligner = config.getString("aligner");
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
	public String aligner;
	/**
	 * Aligner command line argument. String substitution assumes parameters are
	 * fastq, reference, thread count
	 * @return Aligner command line to execute
	 */
	public List<String> getAlignerCommandLine() {
		if (aligner == null) return null;
		switch (aligner) {
			case "bowtie2":
				return ImmutableList.of(
						"bowtie2",
						"--threads",
						"%3$d",
						"--local",
						"--mm",
						"--reorder",
						"-x",
						"%2$s",
						"-U",
						"%1$s");
			case "bwa":
				return ImmutableList.of(
						"bwa",
						"mem",
						"-t",
						"%3$d",
						"-M",
						"%2$s",
						"%1$s");
			default:
				return null;
		}
	}
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
