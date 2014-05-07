package au.edu.wehi.socrates.graph;

import au.edu.wehi.socrates.BreakpointInterval;

/**
 * Weighted trapezoid graph vertex.
 * Start and end positions are both inclusive 
 * @author Daniel Cameron
 *
 */
public class TrapezoidGraphNode {
	public final long start1;
	public final long end1;
	public final long start2;
	public final long end2;
	public final double weight;
	public TrapezoidGraphNode(BreakpointInterval loc) {
		this.start1 = loc.start;
		this.end1 = loc.end;
		this.start2 = loc.start2;
		this.end2 = loc.end2;
		this.weight = (float)loc.qual;
	}
	public TrapezoidGraphNode(TrapezoidGraphNode loc, double additionalWeight) {
		this.start1 = loc.start1;
		this.end1 = loc.end1;
		this.start2 = loc.start2;
		this.end2 = loc.end2;
		this.weight = loc.weight + additionalWeight;
	}
	public TrapezoidGraphNode(long start1, long end1, long start2, long end2, double weight) {
		this.start1 = start1;
		this.end1 = end1;
		this.start2 = start2;
		this.end2 = end2;
		this.weight = weight;
	}
	@Override
	public String toString() {
		return String.format("(%, %, %, %, %)", start1, end1, start2, end2, weight);
	}
}
