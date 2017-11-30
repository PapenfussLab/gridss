package scambler;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Longs;

import au.edu.wehi.idsv.LinearGenomicCoordinate;

public class SgNode {
	public final long inferredPosition;
	public final Read read;
	public SgNode(LinearGenomicCoordinate lgc, Read read, int offset) {
		this.inferredPosition = lgc.getLinearCoordinate(read.getRead().getReferenceIndex(), read.getRead().getUnclippedStart()) + offset;
		this.read = read;
	}
	public State state;
	public TransitiveReductionMark mark;
	public int waitingOnReduction;
	public int waitingOnEmission;
	public enum State {
		UnderConstruction,
		AwaitingReduction,
		Reduced,
	}
	public enum TransitiveReductionMark {
		Vacant,
		InPlay,
		Eliminated,
	}
	/**
	 * Determines whether this node can be compressed into an edge
	 * @return
	 */
	public boolean canCompress() {
		return in.size() != 1 && out.size() == 1;
	}
	public List<SgEdge> in = new ArrayList<>(4);
	public List<SgEdge> out = new ArrayList<>(4);
	//public List<SgEdge> in;
	public static final Ordering<SgNode> ByInferredPosition = new Ordering<SgNode>() {
		@Override
		public int compare(SgNode left, SgNode right) {
			return Longs.compare(left.inferredPosition, right.inferredPosition);
		}
	};
}
