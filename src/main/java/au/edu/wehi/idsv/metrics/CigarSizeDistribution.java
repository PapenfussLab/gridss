package au.edu.wehi.idsv.metrics;

import htsjdk.samtools.CigarOperator;

import java.util.List;

import com.google.api.client.util.Lists;
import com.google.common.collect.Iterables;

/**
 * Distribution of cigar operators
 * @author Daniel Cameron
 *
 */
public class CigarSizeDistribution {
	private final double[][] phred;
	public CigarSizeDistribution(List<CigarDetailMetrics> list) {
		this.phred = new double[CigarOperator.values().length][];
		for (CigarOperator op : CigarOperator.values()) {
			char charop = (char)CigarOperator.enumToCharacter(op);
			this.phred[CigarOperator.enumToBinary(op)] = calcPhred(Lists.newArrayList(Iterables.filter(list, cdm -> cdm.OPERATOR == charop)));
		}
	}
	private static double[] calcPhred(List<CigarDetailMetrics> sc) {
		long total = 0;
		for (CigarDetailMetrics m : sc) {
			total += m.COUNT;
		}
		double[] phred = new double[sc.size()];
		double score = 0;
		long cumsum = total;
		for (int i = 0; i < sc.size(); i++) {
			CigarDetailMetrics m = sc.get(i);
			assert(m.LENGTH == i);
			if (cumsum > 0) {
				score = -10 * Math.log10((double)cumsum / (double)total);
			}
			phred[i] = score;
			cumsum -= m.COUNT;
		}
		if (phred.length == 0) {
			phred = new double[] { 0 };
		}
		return phred;
	}
	/**
	 * Returns the phred scaled probability of a soft clip of at least this length
	 * @param softClipLength
	 * @return soft clip phred score
	 */
	public double getPhred(CigarOperator operator, int length) {
		if (length < 0) return 0;
		double[] opPhred = phred[CigarOperator.enumToBinary(operator)];
		assert(opPhred != null);
		if (length >= opPhred.length) return opPhred[opPhred.length - 1];
		return opPhred[length];
	}
}
