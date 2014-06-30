package au.edu.wehi.idsv;

import java.util.Arrays;

import au.edu.wehi.idsv.vcf.VcfAttributes;

/**
 * Metrics associated with breakpoint evidence
 * @author Daniel Cameron
 *
 */
public class EvidenceMetrics {
	private final int[] v;
	private float baseScore = 0;
	public int get(VcfAttributes evidence) {
		return v[evidence.ordinal()];
	}
	public void set(VcfAttributes evidence, int value) {
		v[evidence.ordinal()] = value;
	}
	@Override
	protected EvidenceMetrics clone() {
		return new EvidenceMetrics(this);
	}
	public EvidenceMetrics(EvidenceMetrics evidence) {
		v = Arrays.copyOf(evidence.v, evidence.v.length);
		baseScore = evidence.baseScore;
	}
	public EvidenceMetrics(float score) {
		this();
		baseScore = score;
	}
	public EvidenceMetrics() {
		 v = new int[VcfAttributes.ASSEMBLY_READS.ordinal() + 1];
	}
	/**
	 * Adds the given additional support to the current evidence
	 * @param evidence evidence to merge with existing evidence 
	 */
	public void add(EvidenceMetrics evidence) {
		if (evidence == null) return;
		for (int i = 0; i < v.length; i++) {
			v[i] += evidence.v[i];
		}
		baseScore += evidence.baseScore;
	}
	/**
	 * Removes the given evidence from the evidence node
	 * @param evidence evidence to remove 
	 */
	public void remove(EvidenceMetrics evidence) {
		if (evidence == null) return;
		for (int i = 0; i < v.length; i++) {
			v[i] -= evidence.v[i];
			//if (v[i] < 0) {
			//	throw new IllegalArgumentException(String.format("Sanity check failure: unable to remove evidence %s from %s.", evidence, this));
			//}
		}
		baseScore -= evidence.baseScore;
	}
	public float getScore() {
		return baseScore
				+ get(VcfAttributes.ASSEMBLY_READS)
				+ get(VcfAttributes.DISCORDANT_READ_PAIR_COUNT)
				+ get(VcfAttributes.UNMAPPED_MATE_READ_COUNT)
				+ get(VcfAttributes.SOFT_CLIP_READ_COUNT);
	}
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder("score=" + Float.toString(getScore()));
		for (VcfAttributes a : VcfAttributes.evidenceValues()) {
			if (v[a.ordinal()] != 0) {
				sb.append(' ');
				sb.append(a.attribute());
				sb.append('=');
				sb.append(v[a.ordinal()]);
			}
		}
		return sb.toString();
	}
}
