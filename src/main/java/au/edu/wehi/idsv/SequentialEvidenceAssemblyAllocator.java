package au.edu.wehi.idsv;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import com.google.common.collect.Iterators;
import com.google.common.collect.Multimap;
import com.google.common.collect.PeekingIterator;
import com.google.common.collect.TreeMultimap;

import htsjdk.samtools.SAMRecord;

public class SequentialEvidenceAssemblyAllocator implements Iterator<SequentialEvidenceAssemblyAllocator.BreakendAssemblyEvidenceSupport> {
	private final LinearGenomicCoordinate lgc;
	private final int maxWindowSize;
	private final PeekingIterator<? extends DirectedEvidence> evidenceIt;
	private final PeekingIterator<SAMRecord> assemblyIt;
	private final ArrayDeque<BreakendAssemblyEvidenceSupport> assemblyBuffer = new ArrayDeque<>();
	private final Multimap<String, BreakendAssemblyEvidenceSupport> evidenceIdToAssembly = TreeMultimap.create();
	public static class BreakendAssemblyEvidenceSupport implements Comparable<BreakendAssemblyEvidenceSupport> {
		public final SAMRecord assemblyRecord;
		public final List<DirectedEvidence> support = new ArrayList<>();
		private BreakendAssemblyEvidenceSupport(SAMRecord assembly) {
			this.assemblyRecord = assembly;
		}
		@Override
		public int compareTo(BreakendAssemblyEvidenceSupport o) {
			return assemblyRecord.getReadName().compareTo(o.assemblyRecord.getReadName());
		}
	}
	/**
	 * Creates an evidence allocator
	 * @param context processing context
	 * @param calls variant calls ordered by position
	 * @param evidence evidence order by breakend position
	 * @param maxCallWindowSize
	 * @param assignEvidenceToSingleBreakpoint uniquely assign evidence to only the highest scoring call
	 */
	public SequentialEvidenceAssemblyAllocator(
			LinearGenomicCoordinate lgc,
			Iterator<? extends DirectedEvidence> reads,
			Iterator<SAMRecord> assemblies,
			int maxAssemblyWindowSize) {
		this.lgc = lgc;
		this.maxWindowSize = maxAssemblyWindowSize;
		this.evidenceIt = Iterators.peekingIterator(reads);
		this.assemblyIt = Iterators.peekingIterator(assemblies);
	}
	@Override
	public boolean hasNext() {
		ensureAllEvidenceAddedToHeadAssembly();
		return !assemblyBuffer.isEmpty();
	}
	@Override
	public BreakendAssemblyEvidenceSupport next() {
		if (!hasNext()) throw new NoSuchElementException();
		BreakendAssemblyEvidenceSupport node = assemblyBuffer.pop();
		for (String evidenceid : new AssemblyAttributes(node.assemblyRecord).getEvidenceIDs()) {
			evidenceIdToAssembly.remove(evidenceid, node);
		}
		return node;
	}
	/**
	 * Invariant: reads that have been iterated over have all possible assemblies loaded
	 */
	private void ensureAllEvidenceAddedToHeadAssembly() {
		// make sure we have an assembly
		if (assemblyBuffer.isEmpty() && assemblyIt.hasNext()) {
			loadAssembly(assemblyIt.next());
		}
		if (assemblyBuffer.isEmpty()) {
			// no more assemblies
			return;
		}
		while (evidenceIt.hasNext() && !isAfterAssembly(assemblyBuffer.peek().assemblyRecord, evidenceIt.peek())) {
			DirectedEvidence read = evidenceIt.next();
			loadAssembliesBefore(lgc.getEndLinearCoordinate(read.getBreakendSummary()) + maxWindowSize);
			allocateRead(read);
		}
	}
	private void allocateRead(DirectedEvidence read) {
		String evidenceid = read.getEvidenceID();
		for (BreakendAssemblyEvidenceSupport node : evidenceIdToAssembly.get(evidenceid)) {
			node.support.add(read);
		}
	}
	private boolean isAfterAssembly(SAMRecord assembly, DirectedEvidence read) {
		return read.getBreakendSummary().referenceIndex > assembly.getReferenceIndex() ||
				(read.getBreakendSummary().referenceIndex == assembly.getReferenceIndex() &&
				read.getBreakendSummary().start - maxWindowSize > assembly.getUnclippedEnd());
	}
	private void loadAssembliesBefore(long position) {
		while (assemblyIt.hasNext() && lgc.getStartLinearCoordinate(assemblyIt.peek()) < position) {
			loadAssembly(assemblyIt.next());
		}
	}
	private void loadAssembly(SAMRecord assembly) {
		BreakendAssemblyEvidenceSupport node = new BreakendAssemblyEvidenceSupport(assembly);
		AssemblyAttributes attr = new AssemblyAttributes(assembly);
		for (String evidenceid : attr.getEvidenceIDs()) {
			evidenceIdToAssembly.put(evidenceid, node);
		}
		assemblyBuffer.add(node);
	}
}
