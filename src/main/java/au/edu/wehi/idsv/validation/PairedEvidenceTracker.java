package au.edu.wehi.idsv.validation;

import au.edu.wehi.idsv.DirectedBreakpoint;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.util.MessageThrottler;
import com.google.common.collect.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.*;

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
	private final HashMap<String, DirectedBreakpoint> unpaired = new HashMap<>();
	private boolean closed = false;
	private HashSet<String> encountered = new HashSet<>();
	private int errors = 0;
	private Writer readNameWriter;

    public PairedEvidenceTracker(String sanitycheck, CloseableIterator<T> iterator, Writer readNameWriter) {
    	this(sanitycheck, iterator);
    	this.readNameWriter = readNameWriter;
    }

    public int errorCount() {
		return errors;
	}

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
		T evidence;
		try {
			if (!it.hasNext()) {
				allMatched();
				return endOfData();
			}
			evidence = it.next();
			isValidNext(evidence);
		} catch (IOException e) {
			log.error(e);
			throw new RuntimeException(e);
		}
		return evidence;
	}

	private void logError(DirectedEvidence evidence, String msgText, String msgType) throws IOException {
		if (!MessageThrottler.Current.shouldSupress(log, msgType)) {
			log.error(msgText);
		}
    	errors++;
    	if (readNameWriter != null) {
			readNameWriter.write(evidence.getUnderlyingSAMRecord().getReadName());
			readNameWriter.write('\n');
		}
	}

	private boolean isValidNext(T evidence) throws IOException {
		ProcessingContext context = evidence.getEvidenceSource().getContext();
		String evidenceId = evidence.getEvidenceID();
		if (evidence instanceof DirectedBreakpoint) {
			DirectedBreakpoint e = (DirectedBreakpoint)evidence;
			String partnerId = e.getRemoteEvidenceID();
			if (encountered.contains(evidenceId)) {
				String msg = String.format("%s: encountered %s multiple times.", name, evidenceId);
				logError(evidence, msg, "duplicate evidence");
				return false;
			}
			if (unpaired.containsKey(partnerId)) {
				DirectedBreakpoint partner = unpaired.remove(partnerId); 
				// Breakpoints must be the same
				if (!partner.getBreakendSummary().remoteBreakpoint().equals(e.getBreakendSummary()) ||
						e.getBreakendSummary().nominal != partner.getBreakendSummary().nominal2 ||
						e.getBreakendSummary().nominal2 != partner.getBreakendSummary().nominal) {
					String msg = String.format("%s: breakpoints %s and %s differ for evidence pair %s %s", name,
							e.getBreakendSummary().toString(context),
							partner.getBreakendSummary().toString(context),
							e.getEvidenceID(),
							partner.getEvidenceID());
					logError(evidence, msg, "asymmetric breakpoint positions");
					return false;
				}
				// Scores must be the same
				if (partner.getBreakpointQual() != e.getBreakpointQual()) {
					String msg = String.format("%s: scores %f and %f differ for evidence pair %s %s", name,
							e.getBreakpointQual(),
							partner.getBreakpointQual(),
							e.getEvidenceID(),
							partner.getEvidenceID());
					logError(evidence, msg, "asymmetric breakpoint quality");
					return false;
				}
				encountered.add(evidenceId);
			} else {
				if (encountered.contains(evidenceId)) {
					String msg = String.format("Duplicate evidence %s encountered for read %s", name,
							e.getBreakpointQual(),
							e.getEvidenceID());
					logError(evidence, msg, "duplicate evidence");
				}
				unpaired.put(evidenceId, e);
			}
		} else {
			if (encountered.contains(evidenceId)) {
				String msg = String.format("%s: encountered %s multiple times.", name, evidenceId);
				logError(evidence, msg, "duplicate evidence");
				return false;
			}
			encountered.add(evidenceId);
		}
		return true;
	}

	private boolean allMatched() throws IOException {
		List<DirectedBreakpoint> list = new ArrayList<>(unpaired.values());
		list.sort(DirectedBreakpoint.ByStartEnd);
		for (DirectedBreakpoint e : list) {
			ProcessingContext context = e.getEvidenceSource().getContext();
			String msg = String.format("%s (%s, %f) unpaired",
					e.getEvidenceID(),
					e.getBreakendSummary().toString(context),
					e.getBreakpointQual());
			logError(e, msg, "unpaired evidence");
		}
		return unpaired.isEmpty();
	}
}
