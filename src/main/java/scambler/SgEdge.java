package scambler;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;
import com.google.common.primitives.Longs;

import au.edu.wehi.idsv.debruijn.PackedSequence;

/**
 * String Graph edge
 * @author Daniel Cameron
 *
 */
public class SgEdge {
	private static final List<SgEdge> NO_EDGES = ImmutableList.of();
	private SgEdge(SgNode from, SgNode to, PackedSequence seq) {
		this.from = from;
		this.to = to;
		this.seq = seq;
		from.out.add(this);
		to.in.add(this);
	}
	/**
	 * Replaces the given compressible node with the equivalent edge
	 * @param node
	 */
	public SgEdge(SgNode node) {
		assert(node.isCompressible());
		SgEdge in = node.in.get(0);
		SgEdge out = node.out.get(0);
		this.from = in.from;
		this.to = out.to;
		this.seq = new PackedSequence(in.seq, out.seq);
		node.edgeState = SgNode.EDGES_FINALISED;
		node.in = NO_EDGES;
		node.out = NO_EDGES;
	}
	/**
	 * Remove the given edge due as it is redundant.
	 */
	public void reduce() {
		from.out.remove(this);
		to.in.remove(this);
	}
	public static List<SgEdge> create(Overlap o) {
		if (o.read2StartRelativeToRead1 < 0) {
			return create(new Overlap(o.read2, o.read1, -o.read2StartRelativeToRead1));
		}
		SgEdge beforeBases = new SgEdge(o.read1.getStartNode(), o.read2.getStartNode(), new PackedSequence(o.read1.getSeq(), 0, o.read1.getSeq().length() - o.overlap));
		SgEdge commonBases = new SgEdge(o.read2.getStartNode(), o.read1.getEndNode(), new PackedSequence(o.read2.getSeq(), 0, o.overlap));
		SgEdge afterBases = new SgEdge(o.read1.getEndNode(), o.read2.getEndNode(), new PackedSequence(o.read2.getSeq(), o.overlap, o.read2.getSeq().length() - o.overlap));
		List<SgEdge> list = new ArrayList<>(3);
		list.add(beforeBases);
		list.add(commonBases);
		list.add(afterBases);
		return list;
	}
	public final SgNode from;
	public final SgNode to;
	public final PackedSequence seq;
	public int length() {
		return seq.length();
	}
	public static final Ordering<SgEdge> BySequenceLength = new Ordering<SgEdge>() {
		@Override
		public int compare(SgEdge left, SgEdge right) {
			return Ints.compare(left.length(), right.length());
		}
	};
	public static final Ordering<SgEdge> ByMaxInferredPosition = new Ordering<SgEdge>() {
		@Override
		public int compare(SgEdge left, SgEdge right) {
			return Longs.compare(
					Math.max(left.from.inferredPosition + left.length(), left.to.inferredPosition),
					Math.max(right.from.inferredPosition + right.length(), right.to.inferredPosition));
		}
	};
	public static final Ordering<SgEdge> ByMinInferredPosition = new Ordering<SgEdge>() {
		@Override
		public int compare(SgEdge left, SgEdge right) {
			return Longs.compare(
					Math.min(left.to.inferredPosition - left.length(), left.from.inferredPosition),
					Math.min(right.to.inferredPosition - left.length(), right.from.inferredPosition));
		}
	};
}
