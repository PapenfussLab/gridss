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
	private DirectedEvidence lastEmitted = null;
	public OrthogonalEvidenceIterator(LinearGenomicCoordinate linear, Iterator<DirectedEvidence> it, int assemblyWindowSize) {
		this.linear = linear;
		this.it = it;
		this.assemblyWindowSize = assemblyWindowSize;
	}
	private boolean outputReady() {
		return !outputQueue.isEmpty() && linear.getStartLinearCoordinate(outputQueue.peek().getBreakendSummary()) < linearPosition - getWindowSize();
	}
	private DirectedEvidence emitNext() {
		DirectedEvidence nextOutput = outputQueue.poll();
		if (lastEmitted != null && DirectedEvidence.ByStartEnd.compare(lastEmitted, nextOutput) > 0) {
			throw new IllegalStateException(String.format("Unable to sort output with window size of %d. %s emitted before %s", assemblyWindowSize, lastEmitted, nextOutput));
		}
		lastEmitted = nextOutput;
		return nextOutput;
	}
	@Override
	protected DirectedEvidence computeNext() {
		if (outputReady()) return emitNext();
		while (it.hasNext()) {
			DirectedEvidence e = it.next();
			add(e);
			call();
			if (outputReady()) return emitNext();
		}
		// reached end of input - flush remaining
		linearPosition = Long.MAX_VALUE;
		call();
		if (outputReady()) return emitNext();
		return endOfData();
	}
	private int getEvidenceWindowSize() { return assemblyWindowSize; }
	private int getAssemblyWindowSize() { return 2 * getEvidenceWindowSize() + 1; }
	private int getWindowSize() { return getAssemblyWindowSize(); }
	/**
	 * Call evidence now outside the merge window
	 */
	private void call() {
		// call evidence that we are guaranteed to have the assembly of in our window
		//  <-----assembly window-----> <-----assembly window----->
		//                             ^
		//                           evidence - check if we're part of any of the assemblies within the window either side 
		//
		// ^ can remove assemblies before this point
		while (!nonAssemblyQueue.isEmpty() && linear.getStartLinearCoordinate(nonAssemblyQueue.peek().getBreakendSummary()) + getEvidenceWindowSize() < linearPosition) {
			outputOrthogonalNonAssembly(nonAssemblyQueue.poll());
		}
		while (!assemblyLookup.isEmpty() && linear.getStartLinearCoordinate(assemblyLookup.peek().getBreakendSummary()) + getAssemblyWindowSize() < linearPosition) {
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
			AssemblyEvidence a = (AssemblyEvidence)e;
			if (!a.isAssemblyFiltered()) {
				assemblyLookup.add(a);
			} else {
				// don't need merge evidence into filtered assemblies 
				outputOrthogonalAssembly(a);
			}
		} else {
			nonAssemblyQueue.add(e);
		}
	}
}
