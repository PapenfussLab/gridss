package au.edu.wehi.socrates.graph;

import java.util.Iterator;

import au.edu.wehi.socrates.BreakpointInterval;
import au.edu.wehi.socrates.BreakpointLocation;
import au.edu.wehi.socrates.DirectedBreakpoint;
import au.edu.wehi.socrates.DirectedEvidence;

import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
/**
 * Calls breakpoints from the given evidence
 * 
 * @author Daniel Cameron
 */
public class EvidenceClusterProcessor implements Iterable<DirectedBreakpoint> {
	// Traverse graph and find all maximum cliques
		// Boxicity = 2
		
		// geometric equivalence:
				// maximal clique == area (ie all points (p1, p2)) in trapezoid of overlap 
				// maximal clique equivalent to all reads matching: a given point
				
				// OEA evidence interval for fusion requires mate to overlap actual fusion sequence [minFragmentSize, maxFragmentSize]
				// OEA evidence interval is [0, maxFragmentSize]
	public EvidenceClusterProcessor() {
	}
	public void addEvidence(Iterator<DirectedEvidence> it) {
		
	}
	@Override
	public Iterator<DirectedBreakpoint> iterator() {
		// TODO Auto-generated method stub
		return null;
	}
}
 