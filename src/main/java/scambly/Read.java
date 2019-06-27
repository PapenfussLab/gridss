package scambly;

import au.edu.wehi.idsv.LinearGenomicCoordinate;
import au.edu.wehi.idsv.debruijn.PackedSequence;
import htsjdk.samtools.SAMRecord;
import org.apache.commons.lang3.tuple.Pair;
import scambler.SgNode;

public class Read {
	private final SAMRecord read;
	public final PackedSequence seq;
	public final long startFirstPosition;
	public final long endFirstPosition;
	public Read(LinearGenomicCoordinate lgc, SAMRecord read, long start, long end) {
		this.read = read;
		this.seq = new PackedSequence(read.getReadBases(), read.getReadNegativeStrandFlag(), read.getReadNegativeStrandFlag());
		this.startFirstPosition = start;
		this.endFirstPosition = end;
	}
	public int readLength() {
		return seq.length();
	}
	@Override
	public String toString() {
		return seq.toString() + " " + read.getReadName();
	}
	public void sanityCheck() {
		for (int i = 0; i < read.getReadLength(); i++) {
			assert(read.getReadBases()[i] == seq.get(i));
		}
		assert(new String(read.getReadBases()).equals(seq.toString()));
	}
}
