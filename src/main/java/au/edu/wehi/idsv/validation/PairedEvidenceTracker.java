package au.edu.wehi.idsv.validation;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

import org.apache.commons.lang3.tuple.Pair;

import com.google.common.collect.AbstractIterator;

import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.DirectedBreakpoint;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.GenomicProcessingContext;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

/**
 * Pass-through iterator that tracks evidence
 * to ensure evidence that should have a remote pairing does so.
 * 
 * @author Daniel Cameron
 *
 */
public class PairedEvidenceTracker<T extends DirectedEvidence> extends AbstractIterator<T> implements CloseableIterator<T> {
	private static final Log log = Log.getInstance(PairedEvidenceTracker.class);
	private final String name;
	private final Iterator<T> it;
	private final GenomicProcessingContext context;
	private final HashMap<String, Pair<BreakpointSummary, Float>> unpaired = new HashMap<>();
	private boolean closed = false;

	public PairedEvidenceTracker(String iteratorName, Iterator<T> it, GenomicProcessingContext context) {
		this.it = it;
		this.name = iteratorName;
		this.context = context;
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
			allMatched();
			return endOfData();
		}
		T evidence = it.next();
		isValidNext(evidence);
		return evidence;
	}

	private boolean isValidNext(T evidence) {
		if (evidence instanceof DirectedBreakpoint) {
			DirectedBreakpoint e = (DirectedBreakpoint)evidence;
			String evidenceId = evidence.getEvidenceID();
			String partnerId = e.getRemoteEvidenceID();
			if (unpaired.containsKey(evidenceId)) {
				String msg = String.format("%s: encountered %s multiple times.", name, evidenceId);
				log.error(msg);
				return false;
			}
			if (unpaired.containsKey(partnerId)) {
				Pair<BreakpointSummary, Float> value = unpaired.remove(partnerId); 
				BreakpointSummary remoteBp = value.getLeft();
				// Breakpoints must be the same
				if (!remoteBp.remoteBreakpoint().equals(e.getBreakendSummary())) {
					String msg = String.format("%s: breakpoints %s and %s differ for evidence pair %s %s", name, e.getBreakendSummary().toString(context), remoteBp.toString(context), e.getEvidenceID(), partnerId); 
					log.error(msg);
					return false;
				}
				// Scores must be the same
				if (value.getRight() != e.getBreakpointQual()) {
					String msg = String.format("%s: scores %f and %f differ for evidence pair %s %s", name, e.getBreakpointQual(), value.getRight(), e.getEvidenceID(), partnerId); 
					log.error(msg);
					return false;
				}
			} else {
				unpaired.put(evidenceId, Pair.of(e.getBreakendSummary(), e.getBreakpointQual()));
			}
		}
		return true;
	}

	private boolean allMatched() {
		for (Entry<String, Pair<BreakpointSummary, Float>> entry : unpaired.entrySet()) {
			log.error(String.format("%s (%s, %f) unpaired",
					entry.getKey(),
					entry.getValue().getLeft().toString(context),
					entry.getValue().getRight()));
		}
		return unpaired.isEmpty();
	}
}
