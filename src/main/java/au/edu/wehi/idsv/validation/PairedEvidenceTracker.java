package au.edu.wehi.idsv.validation;

import java.util.HashMap;
import java.util.Iterator;

import com.google.common.collect.AbstractIterator;

import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.DirectedBreakpoint;
import au.edu.wehi.idsv.DirectedEvidence;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

/**
 * Pass-through iterator that tracks evidence
 * to ensure evidence that should have a remote pairing does so.
 * 
 * Tracking is performed via assertion. If assertions are not enabled,
 * no tracking is performed.
 * 
 * @author Daniel Cameron
 *
 */
public class PairedEvidenceTracker<T extends DirectedEvidence> extends AbstractIterator<T> implements CloseableIterator<T> {
	private static final Log log = Log.getInstance(PairedEvidenceTracker.class);
	/**
	 * Maximum number of unpair records log in full. 
	 */
	private static final int MAX_TO_PRINT = 32;
	private final String name;
	private final Iterator<T> it;
	private boolean closed = false;
	private HashMap<String, BreakpointSummary> unpaired = new HashMap<String, BreakpointSummary>();

	public PairedEvidenceTracker(String iteratorName, Iterator<T> it) {
		this.it = it;
		this.name = iteratorName;
	}

	@Override
	public void close() {
		if (!closed) {
			closed = true;
			CloserUtil.close(it);			
		}
	}

	@Override
	protected T computeNext() {
		if (closed) return endOfData();
		if (!it.hasNext()) {
			assert(allMatched());
			return endOfData();
		}
		T evidence = it.next();
		assert(isValidNext(evidence));
		return evidence;
	}

	private boolean isValidNext(T evidence) {
		if (evidence instanceof DirectedBreakpoint) {
			String evidenceId = evidence.getEvidenceID();
			String partnerId = ((DirectedBreakpoint) evidence).getRemoteEvidenceID();
			if (unpaired.containsKey(evidenceId)) {
				String msg = String.format("%s: encountered %s multiple times.", name, evidenceId);
				log.error(msg);
				return false;
			}
			if (unpaired.containsKey(partnerId)) {
				BreakpointSummary remoteBp = unpaired.remove(partnerId);
				// Breakpoints must be the same
				if (!remoteBp.remoteBreakpoint().equals(((DirectedBreakpoint)evidence).getBreakendSummary())) {
					String msg = String.format("%s: breakpoints %s and %s differ for evidence pair %s %s", name, evidence.getBreakendSummary(), remoteBp, evidence.getEvidenceID(), partnerId); 
					log.error(msg);
					return false;
				}
			} else {
				unpaired.put(evidenceId, (BreakpointSummary)evidence.getBreakendSummary());
			}
		}
		return true;
	}

	private boolean allMatched() {
		if (!unpaired.isEmpty()) {
			String msg = String.format("%s: missing %d evidence pairings: ", name, unpaired.size()); 
			StringBuilder sb = new StringBuilder();
			sb.append(msg);
			int i = MAX_TO_PRINT;
			for (String s : unpaired.keySet()) {
				sb.append(String.format("%s (%s), ", s, unpaired.get(s)));
				if (--i == 0) break;
			}
			if (unpaired.size() > MAX_TO_PRINT) {
				sb.append("...");
			}
			log.error(sb.toString());
		}
		return unpaired.isEmpty();
	}
}
