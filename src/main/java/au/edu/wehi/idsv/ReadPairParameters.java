package au.edu.wehi.idsv;


public class ReadPairParameters {
	/**
	 * Minimum MAPQ of local read anchor to considered as evidence
	 */
	public int minMapq = new SoftClipParameters().minReadMapq;
	public double minAnchorEntropy = new SoftClipParameters().minAnchorEntropy;
	public AdapterHelper adapters = new SoftClipParameters().adapters;
}
