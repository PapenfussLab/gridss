package au.edu.wehi.socrates.graph;

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
	public TrapezoidGraphNode(BreakpointInterval loc) {
		this.startX = loc.start;
		this.endX = loc.end;
		this.startY = loc.start2;
		this.endY = loc.end2;
		this.weight = (float)loc.qual;
	}
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
		return String.format("(x=[%, %], y=[%, %], weight=%)", startX, endX, startY, endY, weight);
	}
}
