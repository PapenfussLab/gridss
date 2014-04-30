package au.edu.wehi.socrates.graph;

import java.util.Iterator;

import au.edu.wehi.socrates.BreakpointInterval;

import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
/**
 * Calculates all maximal cliques from breakpoint interval evidence
 * 
 * @author Daniel Cameron
 */
public class BreakpointIntervalOverlapGraph implements Iterable<BreakpointInterval> {
	private SortedSetMultimap<Long, BreakpointInterval> intervals = TreeMultimap.create(null, null);
	/**
	 * Adds a link between the given intervals
	 * @param interval @see BreakpointInterval evidence to add
	 */
	public void add(BreakpointInterval interval) {
		intervals.put(interval.start1, interval);
	}
	@Override
	public Iterator<BreakpointInterval> iterator() {
		// Traverse graph and find all maximum cliques
		// Boxicity = 2
		
		// geometric equivalence:
				// maximal clique == area (ie all points (p1, p2)) in trapezoid of overlap 
				// maximal clique equivalent to all reads matching: a given point
				
				// OEA evidence interval for fusion requires mate to overlap actual fusion sequence [minFragmentSize, maxFragmentSize]
				// OEA evidence interval is [0, maxFragmentSize]
		return null;
	}
}
 