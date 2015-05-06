package au.edu.wehi.idsv.util;

public class SequenceUtil {
	public static double shannonEntropy(byte[] sequence, int offset, int length) {
		if (length <= 0) return 0;
		if (offset + length > sequence.length) throw new IllegalArgumentException("Sequence too short for specified length and offset");
		int a = 0, c = 0, g = 0, t = 0;
		for (int i = 0; i < length; i++) {
			switch (sequence[offset + i]) {
				case 'A':
				case 'a':
					a++;
					break;
				case 'C':
				case 'c':
					c++;
					break;
				case 'G':
				case 'g':
					g++;
					break;
				case 'T':
				case 't':
					t++;
					break;
			}
		}
		int sum = a + c + g + t;
		double entropy = 0;
		entropy += shannon(a, sum);
		entropy += shannon(c, sum);
		entropy += shannon(g, sum);
		entropy += shannon(t, sum);
		return entropy;
	}
	private static double shannon(int n, int total) {
		if (n == 0) return 0;
		double pr = n / (double)total;
		return -pr * Math.log(pr) / Math.log(2); 
	}
}
