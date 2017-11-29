package scambler;

import java.util.List;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Longs;

import au.edu.wehi.idsv.LinearGenomicCoordinate;

public class SgNode {
	public final long inferredPosition;
	public SgNode(LinearGenomicCoordinate lgc, Read read, int offset) {
		this.inferredPosition = lgc.getLinearCoordinate(read.getRead().getReferenceIndex(), read.getRead().getUnclippedStart()); 
	}
	public List<SgEdge> out;
	//public List<SgEdge> in;
	public static final Ordering<SgNode> ByInferredPosition = new Ordering<SgNode>() {
		@Override
		public int compare(SgNode left, SgNode right) {
			return Longs.compare(left.inferredPosition, right.inferredPosition);
		}
	};
}
