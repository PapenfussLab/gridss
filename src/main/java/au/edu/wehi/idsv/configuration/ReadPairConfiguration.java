package au.edu.wehi.idsv.configuration;

import au.edu.wehi.idsv.AdapterHelper;


public class ReadPairConfiguration {
	/**
	 * Minimum MAPQ of local read anchor to considered as evidence
	 */
	public int minMapq = new SoftClipConfiguration().minReadMapq;
	public double minAnchorEntropy = new SoftClipConfiguration().minAnchorEntropy;
	public AdapterHelper adapters = new SoftClipConfiguration().adapters;
}
