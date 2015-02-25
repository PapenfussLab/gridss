package au.edu.wehi.idsv.sim;

import htsjdk.samtools.reference.ReferenceSequenceFile;

import java.util.NavigableMap;
import java.util.TreeMap;

import com.sun.tools.javac.util.List;

import au.edu.wehi.idsv.DirectedBreakpoint;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;

/**
 * 
 * 
 * @author cameron.d
 *
 */
public class ChromothripticBreakpoints {
	public NavigableMap<Long, DirectedBreakpoint> variants;
	public ChromothripticBreakpoints(ReferenceSequenceFile ref, String chr) {
	}
	/**
	 * Shatter the chromosome
	 */
	public void shatter(ReferenceSequenceFile ref) {
		List<Long> breakpoints;
		// Treat entire genome as one chromosome
		// break at chromosome boundaries
		// shatter chromosomes
		// drop fragments
		// recombine
	}
	public DirectedBreakpoint newRandomVariant(int nUnambiguousPaddingBases, int minDistanceToAdjacentVariant) {
	}
	public int getUnambigousPaddingBaseCount(DirectedBreakpoint breakpoint) {
	}
	public int getUnambigousPaddingBaseCount(DirectedBreakpoint breakpoint) {
	}
}
