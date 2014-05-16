package au.edu.wehi.socrates;

import java.util.Arrays;

import com.sun.tools.javac.code.Attribute.Array;

import au.edu.wehi.socrates.vcf.EvidenceAttributes;

/**
 * Metrics associated with breakpoint evidence
 * @author Daniel Cameron
 *
 */
public class EvidenceMetrics {
	private final int[] v;
	public int get(EvidenceAttributes evidence) {
		return v[evidence.ordinal()];
	}
	public void set(EvidenceAttributes evidence, int value) {
		v[evidence.ordinal()] = value;
	}
	@Override
	protected EvidenceMetrics clone() {
		return new EvidenceMetrics(this);
	}
	public EvidenceMetrics(EvidenceMetrics evidence) {
		v = Arrays.copyOf(evidence.v, evidence.v.length);
	}
	public EvidenceMetrics() {
		 v = new int[EvidenceAttributes.ASSEMBLY_READS.ordinal() + 1];
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
	}
	public double getScore() {
		return get(EvidenceAttributes.ASSEMBLY_READS)
				+ get(EvidenceAttributes.DISCORDANT_READ_PAIR_COUNT)
				+ get(EvidenceAttributes.UNMAPPED_MATE_READ_COUNT)
				+ get(EvidenceAttributes.SOFT_CLIP_READ_COUNT);
	}
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (EvidenceAttributes a : EvidenceAttributes.values()) {
			if (v[a.ordinal()] != 0) {
				sb.append(a.attribute());
				sb.append('=');
				sb.append(v[a.ordinal()]);
				sb.append(' ');
			}
		}
		return sb.toString();
	}
}
