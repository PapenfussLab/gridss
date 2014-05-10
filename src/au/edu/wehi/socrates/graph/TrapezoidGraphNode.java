package au.edu.wehi.socrates.graph;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;

import au.edu.wehi.socrates.BreakpointInterval;

/**
 * Weighted trapezoid graph vertex.
 * Start and end positions are both inclusive 
 * @author Daniel Cameron
 *
 */
public class TrapezoidGraphNode {
	public final long startX;
	public final long endX;
	public final long startY;
	public final long endY;
	public final double weight;
	public TrapezoidGraphNode(TrapezoidGraphNode loc, double additionalWeight) {
		this.startX = loc.startX;
		this.endX = loc.endX;
		this.startY = loc.startY;
		this.endY = loc.endY;
		this.weight = loc.weight + additionalWeight;
	}
	public TrapezoidGraphNode(long startX, long endX, long startY, long endY, double weight) {
		this.startX = startX;
		this.endX = endX;
		this.startY = startY;
		this.endY = endY;
		this.weight = weight;
	}
	@Override
	public String toString() {
		return String.format("(x=[%d, %d], y=[%d, %d], weight=%f)", startX, endX, startY, endY, weight);
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + (int) (endX ^ (endX >>> 32));
		result = prime * result + (int) (endY ^ (endY >>> 32));
		result = prime * result + (int) (startX ^ (startX >>> 32));
		result = prime * result + (int) (startY ^ (startY >>> 32));
		long temp;
		temp = Double.doubleToLongBits(weight);
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
		TrapezoidGraphNode other = (TrapezoidGraphNode) obj;
		if (endX != other.endX)
			return false;
		if (endY != other.endY)
			return false;
		if (startX != other.startX)
			return false;
		if (startY != other.startY)
			return false;
		if (Double.doubleToLongBits(weight) != Double.doubleToLongBits(other.weight))
			return false;
		return true;
	}
	public static Ordering<TrapezoidGraphNode> ByEndXYStartXY = new Ordering<TrapezoidGraphNode>() {
		public int compare(TrapezoidGraphNode o1, TrapezoidGraphNode o2) {
			  return ComparisonChain.start()
			        .compare(o1.endX, o2.endX)
			        .compare(o1.endY, o2.endY)
			        .compare(o1.startX, o2.startX)
			        .compare(o1.startY, o2.startY)
			        .result();
		  }
	};
	public static Ordering<TrapezoidGraphNode> ByStartXYEndXY = new Ordering<TrapezoidGraphNode>() {
		public int compare(TrapezoidGraphNode o1, TrapezoidGraphNode o2) {
			  return ComparisonChain.start()
				  	.compare(o1.startX, o2.startX)
			        .compare(o1.startY, o2.startY)
			        .compare(o1.endX, o2.endX)
			        .compare(o1.endY, o2.endY)
			        .result();
					  
		  }
	};
	public static Ordering<TrapezoidGraphNode> ByEndYXStartYX = new Ordering<TrapezoidGraphNode>() {
		public int compare(TrapezoidGraphNode o1, TrapezoidGraphNode o2) {
			  return ComparisonChain.start()
					 .compare(o1.endY, o2.endY)
					 .compare(o1.endX, o2.endX)
					 .compare(o1.startY, o2.startY)
					 .compare(o1.startX, o2.startX)
					 .result();
		  }
	};
	public static Ordering<TrapezoidGraphNode> ByStartYXEndYX = new Ordering<TrapezoidGraphNode>() {
		public int compare(TrapezoidGraphNode o1, TrapezoidGraphNode o2) {
			return ComparisonChain.start()
					.compare(o1.startY, o2.startY)
					.compare(o1.startX, o2.startX)
					.compare(o1.endY, o2.endY)
					.compare(o1.endX, o2.endX)
					.result();
		  }
	};
}
