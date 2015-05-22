package au.edu.wehi.idsv;

import java.util.Iterator;
import java.util.PriorityQueue;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.HashMultimap;

/**
 * Hydrates assemblies and filters evidence contributing to an assembly
 * @author cameron.d
 *
 */
public class OrthogonalEvidenceIterator extends AbstractIterator<DirectedEvidence> {
	private final Iterator<DirectedEvidence> it;
	private final int assemblyWindowSize;
	private final boolean assembliesOnly;
	private final LinearGenomicCoordinate linear;
	private final PriorityQueue<AssemblyEvidence> assemblyQueue = new PriorityQueue<AssemblyEvidence>(64, DirectedEvidence.ByStartEnd);
	private final PriorityQueue<DirectedEvidence> nonAssemblyQueue = new PriorityQueue<DirectedEvidence>(1024, DirectedEvidence.ByStartEnd); 
	private final PriorityQueue<DirectedEvidence> outputQueue = new PriorityQueue<DirectedEvidence>(1024, DirectedEvidence.ByStartEnd);
	/**
	 * Cache mapping evidence to assembly. Values are either AssemblyEvidence, or List<AssemblyEvidence> depending on how many assemblies
	 * the evidence participates in.
	 */
	private final HashMultimap<String, AssemblyEvidence> assemblyEvidenceLookup = HashMultimap.create(1204, 1);
	private long linearPosition;
	private DirectedEvidence lastEmitted = null;	
	public OrthogonalEvidenceIterator(LinearGenomicCoordinate linear, Iterator<DirectedEvidence> it, int assemblyWindowSize, boolean assembliesOnly) {
		this.linear = linear;
		this.it = it;
		this.assemblyWindowSize = assemblyWindowSize;
		this.assembliesOnly = assembliesOnly;
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
		while (!assemblyQueue.isEmpty() && linear.getStartLinearCoordinate(assemblyQueue.peek().getBreakendSummary()) + getAssemblyWindowSize() < linearPosition) {
			outputOrthogonalAssembly(assemblyQueue.poll());
		}
	}
	private void outputOrthogonalAssembly(AssemblyEvidence e) {
		for (String evidence : e.getEvidenceIDs()) { 
			assemblyEvidenceLookup.remove(evidence, e);
		}
		outputQueue.add(e);
	}
	private void outputOrthogonalNonAssembly(DirectedEvidence e) {
		boolean usedInUnfilteredAssembly = false;
		for (AssemblyEvidence assembly : assemblyEvidenceLookup.get(e.getEvidenceID())) {
			if (assembly instanceof SAMRecordAssemblyEvidence) {
				((SAMRecordAssemblyEvidence)assembly).hydrateEvidenceSet(e);
			}
			usedInUnfilteredAssembly = !assembly.isAssemblyFiltered();
		}
		if (!usedInUnfilteredAssembly && !assembliesOnly) {
			outputQueue.add(e);
		}
	}
	private void add(DirectedEvidence e) {
		linearPosition = linear.getStartLinearCoordinate(e.getBreakendSummary());
		if (e instanceof AssemblyEvidence) {
			AssemblyEvidence a = (AssemblyEvidence)e;
			assemblyQueue.add(a);
			for (String evidence : a.getEvidenceIDs()) {
				assemblyEvidenceLookup.put(evidence, a);
			}
		} else {
			nonAssemblyQueue.add(e);
		}
	}
}
