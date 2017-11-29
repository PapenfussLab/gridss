package scambler;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

import au.edu.wehi.idsv.debruijn.PackedSequence;

/**
 * String Graph edge
 * @author Daniel Cameron
 *
 */
public class SgEdge {
	private SgEdge(SgNode from, SgNode to, PackedSequence seq) {
		this.from = from;
		this.to = to;
		this.seq = seq;
	}
	public static List<SgEdge> create(Overlap o) {
		if (o.read2StartRelativeToRead1 < 0) {
			return create(new Overlap(o.read2, o.read1, -o.read2StartRelativeToRead1));
		}
		SgEdge beforeBases = new SgEdge(o.read1.getStartNode(), o.read2.getStartNode(), new PackedSequence(o.read1.getSeq(), 0, o.read1.getSeq().length() - o.overlap));
		SgEdge commonBases = new SgEdge(o.read2.getStartNode(), o.read1.getEndNode(), new PackedSequence(o.read2.getSeq(), 0, o.overlap));
		SgEdge afterBases = new SgEdge(o.read1.getEndNode(), o.read2.getEndNode(), new PackedSequence(o.read2.getSeq(), o.overlap, o.read2.getSeq().length() - o.overlap));
		o.read1.getStartNode().out.add(beforeBases);
		o.read2.getStartNode().out.add(commonBases);
		o.read1.getEndNode().out.add(afterBases);
		List<SgEdge> list = new ArrayList<>(3);
		list.add(beforeBases);
		list.add(commonBases);
		list.add(afterBases);
		return list;
	}
	SgNode from;
	SgNode to;
	PackedSequence seq;
	public int length() {
		return seq.length();
	}
	public static final Ordering<SgEdge> BySequenceLength = new Ordering<SgEdge>() {
		@Override
		public int compare(SgEdge left, SgEdge right) {
			return Ints.compare(left.seq.length(), right.seq.length());
		}
	};
}
