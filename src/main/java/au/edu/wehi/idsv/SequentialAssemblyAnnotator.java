package au.edu.wehi.idsv;

import java.util.HashMap;
import java.util.Iterator;
import java.util.PriorityQueue;

import com.google.common.collect.AbstractIterator;


/**
 * Rehydrates and annotates assembly evidence 
 * 
 * @author Daniel Cameron
 *
 */
public class SequentialAssemblyAnnotator extends AbstractIterator<SAMRecordAssemblyEvidence> {
	private final LinearGenomicCoordinate lgc;
	private final Iterator<SAMRecordAssemblyEvidence> assit;
	private final Iterator<DirectedEvidence> eit;
	private final int windowSize;
	private final HashMap<String, DirectedEvidence> byEvidenceID = new HashMap<String, DirectedEvidence>();
	private final PriorityQueue<DirectedEvidence> byPosition = new PriorityQueue<DirectedEvidence>(DirectedEvidenceOrder.ByStartEnd);
	private long currentPosition = Integer.MIN_VALUE;
	/**
	 * Creates a new annotator
	 * @param lgc linear coordinate converter
	 * @param assit assembly iterator
	 * @param eit evidence iterator
	 * @param windowSize Size of window to look for supporting evidence.
	 */
	public SequentialAssemblyAnnotator(
			LinearGenomicCoordinate lgc,
			Iterator<SAMRecordAssemblyEvidence> assit,
			Iterator<DirectedEvidence> eit,
			int windowSize) {
		this.lgc = lgc;
		this.assit = assit;
		this.eit = eit;
		this.windowSize = windowSize;
	}
	private void loadUntil(long position) {
		while (eit.hasNext() && currentPosition <= position) {
			DirectedEvidence e = eit.next();
			currentPosition = lgc.getStartLinearCoordinate(e.getBreakendSummary());
			byEvidenceID.put(e.getEvidenceID(), e);
			byPosition.add(e);
		}
	}
	private void flushBefore(long position) {
		while (!byPosition.isEmpty() && lgc.getStartLinearCoordinate(byPosition.peek().getBreakendSummary()) < position) {
			DirectedEvidence e = byPosition.poll();
			byEvidenceID.remove(e.getEvidenceID());
		}
	}
	@Override
	protected SAMRecordAssemblyEvidence computeNext() {
		while (assit.hasNext()) {
			SAMRecordAssemblyEvidence ass = assit.next();
			flushBefore(lgc.getStartLinearCoordinate(ass.getBackingRecord()) - windowSize);
			loadUntil(lgc.getEndLinearCoordinate(ass.getBackingRecord()) + windowSize);
			for(String supportingEvidenceID : ass.getEvidenceIDs()) {
				DirectedEvidence support = byEvidenceID.get(supportingEvidenceID);
				if (support != null) {
					ass.hydrateEvidenceSet(support);
				}
				// evidence should support a single assembly so we're fine to remove it from the lookup once processed
				byEvidenceID.remove(support);
			}
			return ass;
		}
		return endOfData();
	}
}
