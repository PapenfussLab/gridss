package au.edu.wehi.idsv.validation;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.util.MessageThrottler;
import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Ordering;
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
	private final HashMap<String, TrackingInfo> unpaired = new HashMap<>();
	private final boolean trackFullEvidence;
	private final HashSet<String> encountered = new HashSet<>();
	private boolean closed = false;
	private int errors = 0;
	private Writer readNameWriter;
	private ProcessingContext context;

	private static class TrackingInfo {
		public final String evidenceId;
		public final BreakendSummary breakpoint;
		public final float score;
		public String underlyingReadName;
		public final DirectedEvidence fullEvidence;
		public TrackingInfo(DirectedEvidence fullEvidence, boolean trackFullEvidence) {
			this.evidenceId = fullEvidence.getEvidenceID();
			this.breakpoint = fullEvidence.getBreakendSummary();
			this.score = fullEvidence.getBreakendQual();
			this.underlyingReadName = fullEvidence.getUnderlyingSAMRecord().getReadName();
			this.fullEvidence = trackFullEvidence ? fullEvidence : null;
		}
	}

    public PairedEvidenceTracker(String iteratorName, CloseableIterator<T> iterator, Writer readNameWriter, boolean trackFullEvidence) {
    	this(iteratorName, iterator, trackFullEvidence);
    	this.readNameWriter = readNameWriter;
    }

    public int errorCount() {
		return errors;
	}

	public PairedEvidenceTracker(String iteratorName, Iterator<T> it, boolean trackFullEvidence) {
		this.it = it;
		this.name = iteratorName;
		this.trackFullEvidence = trackFullEvidence;
	}

	@Override
	public void close() {
		if (!closed) {
			closed = true;
			unpaired.clear();
			encountered.clear();
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

	private void logError(TrackingInfo ti, String msgText, String msgType) throws IOException {
		if (!MessageThrottler.Current.shouldSupress(log, msgType)) {
			log.error(msgText);
		}
    	errors++;
    	if (readNameWriter != null) {
			readNameWriter.write(ti.underlyingReadName);
			readNameWriter.write('\n');
		}
	}

	private boolean isValidNext(T evidence) throws IOException {
		context = evidence.getEvidenceSource().getContext();
		TrackingInfo ti = new TrackingInfo(evidence, trackFullEvidence);
		String evidenceId = ti.evidenceId;
		if (evidence instanceof DirectedBreakpoint) {
			DirectedBreakpoint e = (DirectedBreakpoint)evidence;
			String partnerId = e.getRemoteEvidenceID();
			if (encountered.contains(evidenceId)) {
				String msg = String.format("%s: encountered %s multiple times.", name, evidenceId);
				logError(ti, msg, "duplicate evidence");
				return false;
			}
			if (unpaired.containsKey(partnerId)) {
				TrackingInfo pti = unpaired.remove(partnerId);
				// Breakpoints must be the same
				BreakpointSummary partnerBreakpoint = (BreakpointSummary)pti.breakpoint;
				if (!partnerBreakpoint.remoteBreakpoint().equals(e.getBreakendSummary()) ||
						e.getBreakendSummary().nominal != partnerBreakpoint.nominal2 ||
						e.getBreakendSummary().nominal2 != partnerBreakpoint.nominal) {
					String msg = String.format("%s: breakpoints %s and %s differ for evidence pair %s %s", name,
							e.getBreakendSummary().toString(context),
							partnerBreakpoint.toString(context),
							e.getEvidenceID(),
							pti.evidenceId);
					logError(ti, msg, "asymmetric breakpoint positions");
					return false;
				}
				// Scores must be the same
				if (pti.score != e.getBreakpointQual()) {
					String msg = String.format("%s: scores %f and %f differ for evidence pair %s %s", name,
							e.getBreakpointQual(),
							pti.score,
							e.getEvidenceID(),
							pti.evidenceId);
					logError(ti, msg, "asymmetric breakpoint quality");
					return false;
				}
				encountered.add(evidenceId);
			} else {
				if (encountered.contains(evidenceId)) {
					String msg = String.format("Duplicate evidence %s encountered for read %s", name,
							e.getBreakpointQual(),
							e.getEvidenceID());
					logError(ti, msg, "duplicate evidence");
				}
				unpaired.put(evidenceId, ti);
			}
		} else {
			if (encountered.contains(evidenceId)) {
				String msg = String.format("%s: encountered %s multiple times.", name, evidenceId);
				logError(ti, msg, "duplicate evidence");
				return false;
			}
			encountered.add(evidenceId);
		}
		return true;
	}

	private boolean allMatched() throws IOException {
		List<TrackingInfo> list = new ArrayList<>(unpaired.values());
		list.sort(BreakendSummary.ByStartEnd.onResultOf(ti -> ti.breakpoint));
		for (TrackingInfo ti : list) {
			String msg = String.format("%s: %s (%s, %f) unpaired", name,
					ti.evidenceId,
					ti.breakpoint.toString(context),
					ti.score);
			logError(ti, msg, "unpaired evidence");
		}
		return unpaired.isEmpty();
	}
}
