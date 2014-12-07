package au.edu.wehi.idsv;

import java.util.Iterator;
import java.util.PriorityQueue;

import com.google.common.collect.AbstractIterator;

/**
 * Filters sorted evidence such that all raw evidence contributing to an assembly is replaced by
 * the assembly 
 * @author cameron.d
 *
 */
public class OrthogonalEvidenceIterator extends AbstractIterator<DirectedEvidence> {
	private final Iterator<DirectedEvidence> it;
	private final int assemblyWindowSize;
	private final LinearGenomicCoordinate linear;
	private final PriorityQueue<AssemblyEvidence> assemblyLookup = new PriorityQueue<AssemblyEvidence>(64, DirectedEvidence.ByStartEnd);
	private final PriorityQueue<DirectedEvidence> nonAssemblyQueue = new PriorityQueue<DirectedEvidence>(1024, DirectedEvidence.ByStartEnd); 
	private final PriorityQueue<DirectedEvidence> outputQueue = new PriorityQueue<DirectedEvidence>(1024, DirectedEvidence.ByStartEnd);
	private long linearPosition;
	public OrthogonalEvidenceIterator(LinearGenomicCoordinate linear, Iterator<DirectedEvidence> it, int assemblyWindowSize) {
		this.linear = linear;
		this.it = it;
		this.assemblyWindowSize = assemblyWindowSize;
	}
	private boolean outputReady() {
		return !outputQueue.isEmpty() && linear.getStartLinearCoordinate(outputQueue.peek().getBreakendSummary()) < linearPosition - assemblyWindowSize;
	}
	@Override
	protected DirectedEvidence computeNext() {
		if (outputReady()) return outputQueue.poll();
		while (it.hasNext()) {
			DirectedEvidence e = it.next();
			add(e);
			call();
			if (outputReady()) return outputQueue.poll();
		}
		// reached end of input - flush remaining
		linearPosition = Long.MAX_VALUE;
		call();
		if (outputReady()) return outputQueue.poll();
		return endOfData();
	}
	/**
	 * Call evidence now outside the merge window
	 */
	private void call() {
		// call evidence that we are guaranteed to have the assembly of in our window
		if (!nonAssemblyQueue.isEmpty() && linear.getStartLinearCoordinate(nonAssemblyQueue.peek().getBreakendSummary()) + assemblyWindowSize < linearPosition) {
			outputOrthogonalNonAssembly(nonAssemblyQueue.poll());
		}
		if (!assemblyLookup.isEmpty() && linear.getStartLinearCoordinate(assemblyLookup.peek().getBreakendSummary()) + 2 * assemblyWindowSize < linearPosition) {
			outputOrthogonalAssembly(assemblyLookup.poll());
		}
	}
	private void outputOrthogonalAssembly(AssemblyEvidence e) {
		outputQueue.add(e);
	}
	private void outputOrthogonalNonAssembly(DirectedEvidence e) {
		for (AssemblyEvidence assembly : assemblyLookup) {
			if (assembly.isPartOfAssemblyBreakend(e)) {
				// evidence is part of our assembly: no need to output
				return;
			}
		}
		outputQueue.add(e);
	}
	private void add(DirectedEvidence e) {
		linearPosition = linear.getStartLinearCoordinate(e.getBreakendSummary());
		if (e instanceof AssemblyEvidence) {
			assemblyLookup.add((AssemblyEvidence)e);
		} else {
			nonAssemblyQueue.add(e);
		}
	}
}
