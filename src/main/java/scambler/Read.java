package scambler;

import org.apache.commons.lang3.tuple.Pair;

import au.edu.wehi.idsv.LinearGenomicCoordinate;
import au.edu.wehi.idsv.debruijn.PackedSequence;
import htsjdk.samtools.SAMRecord;

public class Read {
	private final SAMRecord read;
	private final PackedSequence seq;
	private final SgNode startNode;
	private final SgNode endNode;
	private Read mate;
	private Read(LinearGenomicCoordinate lgc,SAMRecord read) {
		this.read = read;
		this.seq = new PackedSequence(read.getReadBases(), false, false);
		if (read.getReadUnmappedFlag()) {
			throw new IllegalStateException("Orientation of unmapped read unknown");
		}
		this.startNode = new SgNode(lgc, this, 0);
		this.endNode = new SgNode(lgc, this, seq.length());
	}
	public static Read create(LinearGenomicCoordinate lgc, SAMRecord read) {
		return new Read(lgc, read);
	}
	public static Pair<Read, Read> create(LinearGenomicCoordinate lgc,SAMRecord read1, SAMRecord read2) {
		Read rn1 = new Read(lgc, read1);
		Read rn2 = new Read(lgc, read2);
		rn1.mate = rn2;
		rn2.mate = rn1;
		return Pair.of(rn1, rn2);
	}
	public PackedSequence getSeq() {
		return seq;
	}
	public SAMRecord getRead() {
		return read;
	}
	public Read getMate() {
		return mate;
	}
	public SgNode getStartNode() {
		return startNode;
	}
	public SgNode getEndNode() {
		return endNode;
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
