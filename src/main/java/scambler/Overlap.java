package scambler;

import au.edu.wehi.idsv.debruijn.PackedSequence;

public class Overlap {
	public Overlap(Read read1, Read read2, int read2StartRelativeToRead1) {
		this.read1 = read1;
		this.read2 = read2;
		this.read2StartRelativeToRead1 = read2StartRelativeToRead1;
		this.matchingBases = PackedSequence.overlapMatches(read1.getSeq(), read2.getSeq(), read2StartRelativeToRead1);
		this.overlap = PackedSequence.overlapLength(read1.getSeq(), read2.getSeq(), read2StartRelativeToRead1);
	}
	public final Read read1;
	public final Read read2;
	public final int read2StartRelativeToRead1;
	public final int matchingBases;
	public final int overlap;
	/**
	 * Distance (in base pairs) between the actual or expected reference genome alignment positions
	 * and that inferred by this overlap. 
	 */
	public int deviationFromAlignment() {
		// TODO: calculate
		return 0;
	}
}
