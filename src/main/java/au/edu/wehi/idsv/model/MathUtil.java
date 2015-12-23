package au.edu.wehi.idsv.model;

public class MathUtil {
	public static double phredToPr(double phred) {
		return Math.pow(10, -phred / 10);
	}
	public static double prToPhred(double pr) {
		if (pr == 0) return 0;
		return -10 * Math.log10(pr);
	}
	public static double mapqToPr(double mapq) {
		return 1 - MathUtil.phredToPr(mapq);
	}
	public static double mapqToPr(int mapq1, int mapq2) {
		return mapqToPr(mapq1) * mapqToPr(mapq2);
	}
}
