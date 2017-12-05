package scambler;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Longs;

import au.edu.wehi.idsv.LinearGenomicCoordinate;

public class SgNode {
	public List<SgEdge> in = new ArrayList<>(4);
	public List<SgEdge> out = new ArrayList<>(4);
	public final long inferredPosition;
	public final Read read;
	public State state = State.UnderConstruction;
	public TransitiveReductionMark mark = TransitiveReductionMark.Vacant; 
	public int waitingOnReduction;
	public int waitingOnEmission;
	public enum State {
		UnderConstruction,
		AwaitingReduction,
		Reduced,
		EmittedFromReduction,
		/**
		 * This node has been removed from the graph
		 */
		Removed,
	}
	public enum TransitiveReductionMark {
		Vacant,
		InPlay,
		Eliminated,
	}
	public SgNode(LinearGenomicCoordinate lgc, Read read, int offset) {
		this.inferredPosition = lgc.getLinearCoordinate(read.getRead().getReferenceIndex(), read.getRead().getUnclippedStart()) + offset;
		this.read = read;
	}
	/**
	 * Determines whether this node can be compressed into an edge
	 * @return
	 */
	public boolean canCompress() {
		return in.size() != 1 && out.size() == 1;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder(String.format("\"%s\"(%d){", read.getRead().getReadName(), inferredPosition));
		for (SgEdge e : out) {
			sb.append(String.format("%s->\"%s\"(%d),", e.seq, e.to.read.getRead().getReadName(), e.to.inferredPosition));
		}
		sb.append("}");
		return sb.toString();
	}
	//public List<SgEdge> in;
	public static final Ordering<SgNode> ByInferredPosition = new Ordering<SgNode>() {
		@Override
		public int compare(SgNode left, SgNode right) {
			return Longs.compare(left.inferredPosition, right.inferredPosition);
		}
	};
}
