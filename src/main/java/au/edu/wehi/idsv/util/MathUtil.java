package au.edu.wehi.idsv.util;

public class MathUtil {
	public static double phredToPr(double phred) {
		return StrictMath.pow(10, -phred / 10);
	}
	public static double prToPhred(double pr) {
		if (pr == 0) return 0;
		return -10 * StrictMath.log10(pr);
	}
	/**
	 * Returns the phred-scaled error probability of no
	 * errors given a set of independent error sources
	 * @param score phred-scaled error probabilities of each error source
	 * @return
	 */
	public static double phredOr(double... phredScore) {
		assert(phredScore.length > 0);
		double resultPhred = phredScore[0];
		// -10 * log10(1 - (1 - 10^(-a/10)) * (1 - 10^(-b/10)))
		// = -10 * log10(10^(-a/10)+10^(-b/10)-10^((-a-b)/10))
		for (int i = 1; i < phredScore.length; i++) {
			resultPhred = prToPhred(phredToPr(resultPhred) + phredToPr(phredScore[i]) - phredToPr(resultPhred + phredScore[i]));
		}
		return resultPhred;
	}
	/**
	 * Computes the average of the given integers without being vulnerable to integer overflow
	 * @param values
	 * @return
	 */
	public static int average(int... values) {
		long sum = 0;
		for (int v : values) {
			sum += v;
		}
		return (int)(sum / values.length);
	}
}
