package au.edu.wehi.idsv.configuration;

import org.apache.commons.configuration.Configuration;

public class AnchorRealignmentConfiguration {
	public static final String CONFIGURATION_PREFIX = "anchorRealignment";
	public AnchorRealignmentConfiguration(Configuration contig) {
		contig = contig.subset(CONFIGURATION_PREFIX);
		perform = contig.getBoolean("perform");
		realignmentMinimumAnchorRetainment = contig.getFloat("realignmentWindowReadLengthMultiples");
		realignmentWindowReadLengthMultiples = contig.getDouble("realignmentMinimumAnchorRetainment");
	}
	/**
	 * Perform local Smith-Waterman realignment of assemblies
	 */
	public boolean perform = true;
	/**
	 * Local realignment window size (in multiple of read length).
	 * This window is in addition to the 
	 */
	public double realignmentWindowReadLengthMultiples = 0.2;
	/**
	 * Minimum portion of the anchored bases that remained anchored
	 * 
	 * When performing realignment around small indels and tandem duplications,
	 * if the breakend sequence is significantly longer than the anchor sequence
	 * then realignment will align the breakend to the other side of the SV, and
	 * consider the anchor to be the new breakend sequence. To prevent this,
	 * realignment is aborted if more than this portion of the initial reference
	 * bases are no longer aligned to the reference after realignment.
	 * 
	 */
	public float realignmentMinimumAnchorRetainment = 0.5f;
}