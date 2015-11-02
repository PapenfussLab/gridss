package au.edu.wehi.idsv.sim;

/**
 * Places variants sequentially, with a padding distance between them.
 * @author Daniel Cameron
 *
 */
public class SequentialVariantPlacer {
	@SuppressWarnings("serial")
	public static class ContigExhaustedException extends Exception {
		public ContigExhaustedException(String msg) {
			super(msg);
		}
	}
	private int position = 0;
	private int distanceBetweenVariants = 2000; // 2 * 500bp fragments buffer on either size
	private int minimumDistanceFromN = -1;
	private byte[] reference;
	public SequentialVariantPlacer(byte[] reference) {
		this.reference = reference;
	}
	public SequentialVariantPlacer(byte[] reference, int distanceBetweenVariants) {
		this.reference = reference;
		this.distanceBetweenVariants = distanceBetweenVariants;
		this.minimumDistanceFromN = distanceBetweenVariants;
	}
	public SequentialVariantPlacer(byte[] reference, int distanceBetweenVariants, int minimumDistanceFromN) {
		this.reference = reference;
		this.distanceBetweenVariants = distanceBetweenVariants;
		this.minimumDistanceFromN = minimumDistanceFromN;
	}
	public int getNext(int featureSize) throws ContigExhaustedException {
		if (featureSize <= 0) throw new IllegalArgumentException("Feature size must be at least 1 base");
		int minN = getMinimumDistanceFromN();
		
		// find the next suitable location to place the given feature
		position = position + getDistanceBetweenVariants() + 1;
		if (getMinimumDistanceFromN() != 0) {
			int offset = getOffsetOfNextUnambiguousSequence(Math.max(0, position - 1 - minN), minN * 2 + featureSize);
			if (offset < 0) {
				throw new ContigExhaustedException("Unable to get next variant location: no unambigious locations remains (contig exhausted).");
			}
			position = offset + 1 + minN;
		}
		if (position + featureSize + minN > reference.length) {
			throw new ContigExhaustedException("Unable to get next variant location: no valid location remains (contig exhausted).");
		}
		int startPosition = position;
		position = position + featureSize - 1;
		return startPosition;
	}
	private int getOffsetOfNextUnambiguousSequence(int offset, int length) {
		for (int i = offset; i < reference.length && i < offset + length; i++) {
			if (isAmbiguous(reference[i])) {
				offset = i + 1;
			}
		}
		return offset;
	}
	private static boolean isAmbiguous(byte b) {
		return b == 'N';
	}
	/**
	 * Distance between successive variants
	 * @return distance between successive variants
	 */
	public int getDistanceBetweenVariants() {
		return distanceBetweenVariants;
	}
	/**
	 * Distance between successive variants
	 * @return
	 */
	public void setDistanceBetweenVariants(int distanceBetweenVariants) {
		this.distanceBetweenVariants = distanceBetweenVariants;
	}
	/**
	 * Gets the minimum distance a variant can be from an N base in the reference
	 * @return  minimum distance from an N base
	 */
	public int getMinimumDistanceFromN() {
		return minimumDistanceFromN >= 0 ? minimumDistanceFromN : distanceBetweenVariants;
	}
	/**
	 * set the minimum distance a variant can be from an N base in the reference
	 * @param minimumDistanceFromN minimum distance from an N base
	 */
	public void setMinimumDistanceFromN(int minimumDistanceFromN) {
		this.minimumDistanceFromN = minimumDistanceFromN;
	}
}
