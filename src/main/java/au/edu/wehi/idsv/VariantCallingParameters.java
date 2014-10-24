package au.edu.wehi.idsv;

public class VariantCallingParameters {
	/**
	 * Maximum somatic p-value for flagging a variant as somatic
	 */
	public double somaticPvalueThreshold = 0.001;
	/**
	 * Call breakends only on assembled contigs
	 */
	public boolean callOnlyAssemblies = false;
	/**
	 * Minimum indel size
	 */
	public int minIndelSize = 16;
}
