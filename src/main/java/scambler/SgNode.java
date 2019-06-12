package scambler;

import au.edu.wehi.idsv.LinearGenomicCoordinate;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Longs;

import java.util.ArrayList;
import java.util.List;

public class SgNode {
	/**
	 * Edges are still being added to this node
	 */
	public static final int EDGES_UNDER_CONSTRUCTION = 0;
	/**
	 * All edges have been added
	 */
	public static final int EDGES_CONSTRUCTED = 1;
	/**
	 * Node has had transitive reduction (and compression) performed
	 */
	public static final int EDGES_TRANSISTIVE_REDUCED = 2;
	/**
	 * Transitive reduction and node compression has been completed.
	 * No further changes to this node will be made
	 */
	public static final int EDGES_FINALISED = 3;
	public List<SgEdge> in = new ArrayList<>(4);
	public List<SgEdge> out = new ArrayList<>(4);
	public final long inferredPosition;
	public final Read read;
	//public int waitingOnReduction;
	//public int waitingOnEmission;
	public int edgeState = 0;
	public TransitiveReductionMark mark = TransitiveReductionMark.Vacant;
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
	public boolean isCompressible() {
		return in.size() == 1 && out.size() == 1;
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
