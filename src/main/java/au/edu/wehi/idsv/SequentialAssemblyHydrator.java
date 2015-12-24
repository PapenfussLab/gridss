package au.edu.wehi.idsv;

import java.util.ArrayDeque;
import java.util.HashMap;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.google.common.collect.PeekingIterator;


/**
 * Rehydrates assemblies with their supporting evidence 
 * 
 * @author Daniel Cameron
 *
 */
public class SequentialAssemblyHydrator implements PeekingIterator<SAMRecordAssemblyEvidence> {
	private final LinearGenomicCoordinate lgc;
	private final PeekingIterator<SAMRecordAssemblyEvidence> assit;
	private final Iterator<DirectedEvidence> eit;
	private final int windowSize;
	private final HashMap<String, SAMRecordAssemblyEvidence> byEvidenceID = new HashMap<String, SAMRecordAssemblyEvidence>();
	private final Queue<SAMRecordAssemblyEvidence> byPosition = new ArrayDeque<SAMRecordAssemblyEvidence>();
	private long currentPosition = Long.MIN_VALUE;
	/**
	 * Creates a new annotator
	 * @param lgc linear coordinate converter
	 * @param assit assembly iterator
	 * @param eit evidence iterator
	 * @param windowSize Size of window to look for supporting evidence.
	 */
	public SequentialAssemblyHydrator(
			LinearGenomicCoordinate lgc,
			PeekingIterator<SAMRecordAssemblyEvidence> assit,
			Iterator<DirectedEvidence> eit,
			int windowSize) {
		this.lgc = lgc;
		this.assit = assit;
		this.eit = eit;
		this.windowSize = windowSize;
	}
	private void loadAssembliesAround(long position) {
		// since all assemblies are generated from evidence,
		// there isn't a memory consumption issue jumping to
		// the next possible supporting read.
		while (assit.hasNext() && lgc.getStartLinearCoordinate(assit.peek().getBreakendSummary()) <= position + windowSize) {
			SAMRecordAssemblyEvidence ass = assit.next();
			for (String eid : ass.getEvidenceIDs()) {
				byEvidenceID.put(eid, ass);
			}
			byPosition.add(ass);
		}
	}
	private boolean assemblyFullAnnotated() {
		return !byPosition.isEmpty() && lgc.getStartLinearCoordinate(byPosition.peek().getBreakendSummary()) < currentPosition - windowSize;
	}
	private void ensureAnnotated() {
		while (!assemblyFullAnnotated() && eit.hasNext()) {
			DirectedEvidence e = eit.next();
			currentPosition = lgc.getStartLinearCoordinate(e.getBreakendSummary());
			String eid = e.getEvidenceID();
			loadAssembliesAround(currentPosition);
			SAMRecordAssemblyEvidence ass = byEvidenceID.get(eid);
			if (ass != null) {
				ass.hydrateEvidenceSet(e);
				byEvidenceID.remove(eid);
			}
		}
		// No more supporting evidence: just pass through the remaining assemblies
		if (!eit.hasNext() && byPosition.isEmpty() && assit.hasNext()) {
			currentPosition = Long.MAX_VALUE;
			byPosition.add(assit.next());
		}
	}
	@Override
	public boolean hasNext() {
		ensureAnnotated();
		return !byPosition.isEmpty();
	}
	@Override
	public SAMRecordAssemblyEvidence next() {
		ensureAnnotated();
		SAMRecordAssemblyEvidence ass = byPosition.poll();
		for (String eid : ass.getEvidenceIDs()) {
			byEvidenceID.remove(eid);
		}
		return ass;
	}
	@Override
	public SAMRecordAssemblyEvidence peek() {
		ensureAnnotated();
		if (byPosition.isEmpty()) throw new NoSuchElementException();
		return byPosition.peek();
	}
	@Override
	public void remove() {
		throw new UnsupportedOperationException("remove");
	}
}
