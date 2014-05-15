package au.edu.wehi.socrates.graph;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;

import au.edu.wehi.socrates.BreakpointSummary;
import au.edu.wehi.socrates.EvidenceMetrics;

/**
 * Weighted trapezoid graph vertex.
 * Start and end positions are both inclusive 
 * @author Daniel Cameron
 *
 */
public class GraphNode {
	public final long startX;
	public final long endX;
	public final long startY;
	public final long endY;
	public final EvidenceMetrics evidence;
	public GraphNode(GraphNode loc, EvidenceMetrics additionalEvidence) {
		this.startX = loc.startX;
		this.endX = loc.endX;
		this.startY = loc.startY;
		this.endY = loc.endY;
		this.evidence = loc.evidence;
		if (loc.evidence != null) {
			this.evidence.add(additionalEvidence);
		} else {
			this.evidence = additionalEvidence;
		}
	}
	public GraphNode(long startX, long endX, long startY, long endY, EvidenceMetrics evidence) {
		this.startX = startX;
		this.endX = endX;
		this.startY = startY;
		this.endY = endY;
		this.evidence = evidence;
	}
	@Override
	public String toString() {
		return String.format("(x=[%d, %d], y=[%d, %d], %s)", startX, endX, startY, endY, evidence);
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + (int) (endX ^ (endX >>> 32));
		result = prime * result + (int) (endY ^ (endY >>> 32));
		result = prime * result + (int) (startX ^ (startX >>> 32));
		result = prime * result + (int) (startY ^ (startY >>> 32));
		int temp;
		temp = evidence.hashCode();
		result = prime * result + (int) (temp ^ (temp >>> 32));
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		GraphNode other = (GraphNode) obj;
		if (endX != other.endX)
			return false;
		if (endY != other.endY)
			return false;
		if (startX != other.startX)
			return false;
		if (startY != other.startY)
			return false;
		if (!evidence.equals(other.evidence))
			return false;
		return true;
	}
	public static Ordering<GraphNode> ByEndXYStartXY = new Ordering<GraphNode>() {
		public int compare(GraphNode o1, GraphNode o2) {
			  return ComparisonChain.start()
			        .compare(o1.endX, o2.endX)
			        .compare(o1.endY, o2.endY)
			        .compare(o1.startX, o2.startX)
			        .compare(o1.startY, o2.startY)
			        .result();
		  }
	};
	public static Ordering<GraphNode> ByStartXYEndXY = new Ordering<GraphNode>() {
		public int compare(GraphNode o1, GraphNode o2) {
			  return ComparisonChain.start()
				  	.compare(o1.startX, o2.startX)
			        .compare(o1.startY, o2.startY)
			        .compare(o1.endX, o2.endX)
			        .compare(o1.endY, o2.endY)
			        .result();
					  
		  }
	};
	public static Ordering<GraphNode> ByEndYXStartYX = new Ordering<GraphNode>() {
		public int compare(GraphNode o1, GraphNode o2) {
			  return ComparisonChain.start()
					 .compare(o1.endY, o2.endY)
					 .compare(o1.endX, o2.endX)
					 .compare(o1.startY, o2.startY)
					 .compare(o1.startX, o2.startX)
					 .result();
		  }
	};
	public static Ordering<GraphNode> ByStartYXEndYX = new Ordering<GraphNode>() {
		public int compare(GraphNode o1, GraphNode o2) {
			return ComparisonChain.start()
					.compare(o1.startY, o2.startY)
					.compare(o1.startX, o2.startX)
					.compare(o1.endY, o2.endY)
					.compare(o1.endX, o2.endX)
					.result();
		  }
	};
}
