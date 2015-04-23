package au.edu.wehi.idsv.graph;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;

/**
 * Weighted rectangle graph vertex.
 * Start and end positions are both inclusive 
 * @author Daniel Cameron
 *
 */
public class RectangleGraphNode {
	public final long startX;
	public final long endX;
	public final long startY;
	public final long endY;
	public long weight;
	public RectangleGraphNode(long startX, long endX, long startY, long endY, long weight) {
		assert(weight > 0);
		this.startX = startX;
		this.endX = endX;
		this.startY = startY;
		this.endY = endY;
		this.weight = weight;
	}
	public boolean isSameCoordinate(RectangleGraphNode node) {
		return startX == node.startX &&
				endX == node.endX &&
				startY == node.startY &&
				endY == node.endY;
	}
	/**
	 * Flips the axes for the given graph node
	 * @return
	 */
	public RectangleGraphNode flipAxis() {
		return new RectangleGraphNode(startY, endY, startX, endX, weight);
	}
	@Override
	public String toString() {
		return String.format("(x=[%d, %d], y=[%d, %d], %d)", startX, endX, startY, endY, weight);
	}
	public static Ordering<RectangleGraphNode> ByEndXYStartXY = new Ordering<RectangleGraphNode>() {
		public int compare(RectangleGraphNode o1, RectangleGraphNode o2) {
			  return ComparisonChain.start()
			        .compare(o1.endX, o2.endX)
			        .compare(o1.endY, o2.endY)
			        .compare(o1.startX, o2.startX)
			        .compare(o1.startY, o2.startY)
			        .result();
		  }
	};
	public static Ordering<RectangleGraphNode> ByStartXYEndXY = new Ordering<RectangleGraphNode>() {
		public int compare(RectangleGraphNode o1, RectangleGraphNode o2) {
			  return ComparisonChain.start()
				  	.compare(o1.startX, o2.startX)
			        .compare(o1.startY, o2.startY)
			        .compare(o1.endX, o2.endX)
			        .compare(o1.endY, o2.endY)
			        .result();
		  }
	};
	public static Ordering<RectangleGraphNode> ByStartXY = new Ordering<RectangleGraphNode>() {
		public int compare(RectangleGraphNode o1, RectangleGraphNode o2) {
			  return ComparisonChain.start()
				  	.compare(o1.startX, o2.startX)
			        .compare(o1.startY, o2.startY)
			        .result();
		  }
	};
	public static Ordering<RectangleGraphNode> ByEndYXStartYX = new Ordering<RectangleGraphNode>() {
		public int compare(RectangleGraphNode o1, RectangleGraphNode o2) {
			  return ComparisonChain.start()
					 .compare(o1.endY, o2.endY)
					 .compare(o1.endX, o2.endX)
					 .compare(o1.startY, o2.startY)
					 .compare(o1.startX, o2.startX)
					 .result();
		  }
	};
	public static Ordering<RectangleGraphNode> ByStartYXEndYX = new Ordering<RectangleGraphNode>() {
		public int compare(RectangleGraphNode o1, RectangleGraphNode o2) {
			return ComparisonChain.start()
					.compare(o1.startY, o2.startY)
					.compare(o1.startX, o2.startX)
					.compare(o1.endY, o2.endY)
					.compare(o1.endX, o2.endX)
					.result();
		  }
	};
	public static Ordering<RectangleGraphNode> ByEndXStartYEndY = new Ordering<RectangleGraphNode>() {
		public int compare(RectangleGraphNode o1, RectangleGraphNode o2) {
			return ComparisonChain.start()
					.compare(o1.endX, o2.endX)
					.compare(o1.startY, o2.startY)
					.compare(o1.endY, o2.endY)
					.result();
		  }
	};
	public static Ordering<RectangleGraphNode> ByEndY = new Ordering<RectangleGraphNode>() {
		public int compare(RectangleGraphNode o1, RectangleGraphNode o2) {
			return ComparisonChain.start()
					.compare(o1.endY, o2.endY)
					.result();
		  }
	};
}
