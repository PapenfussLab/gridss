package au.edu.wehi.idsv.sim;

import java.util.Random;

public class RandomBaseGenerator implements RandomSequenceGenerator {
	private final Random random;
	private static final byte[] DNA_BASES = {'A', 'T', 'C', 'G', };
	public RandomBaseGenerator(int seed) {
		random = new Random(seed);
	}
	public byte[] getBases(int length) {
		byte[] bases = new byte[length];
		for (int i = 0; i < length; i++) {
			bases[i] = DNA_BASES[random.nextInt(4)];
		}
		return bases;
	}
}
