package au.edu.wehi.idsv;

import java.util.Iterator;
import java.util.Set;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;
import com.google.common.collect.Sets;

/**
 * Filters directed evidence to only contain evidence involving the
 * given pair of chromosomes
 * @author Daniel Cameron
 *
 */
public class DirectedEvidenceChromosomePairFilter extends AbstractIterator<DirectedEvidence> {
	private final PeekingIterator<DirectedEvidence> it;
	private final Set<Integer> allowed = Sets.newHashSet();
	public DirectedEvidenceChromosomePairFilter(Iterator<DirectedEvidence> underlying, int referenceIndex1, int referenceIndex2) {
		this.it = Iterators.peekingIterator(underlying);
		this.allowed.add(referenceIndex1);
		this.allowed.add(referenceIndex2);
	}
	private boolean filterOut(DirectedEvidence evidence) {
		BreakendSummary loc = evidence.getBreakendSummary();
		if (allowed.contains(loc.referenceIndex)) return false;
		if (loc instanceof BreakpointSummary && allowed.contains(((BreakpointSummary)loc).referenceIndex2)) return false;
		return true;
	}
	@Override
	protected DirectedEvidence computeNext() {
		while (it.hasNext() && filterOut(it.peek())) it.next();
		if (!it.hasNext()) return endOfData();
		return it.next();
	} 
}
