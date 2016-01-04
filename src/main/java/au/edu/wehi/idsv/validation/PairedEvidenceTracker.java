package au.edu.wehi.idsv.validation;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import java.util.HashMap;
import java.util.Iterator;

import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.DirectedBreakpoint;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.DiscordantReadPair;
import au.edu.wehi.idsv.MaximalEvidenceCliqueIterator;
import au.edu.wehi.idsv.RealignedRemoteSAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.RealignedRemoteSoftClipEvidence;
import au.edu.wehi.idsv.RealignedSAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.RealignedSoftClipEvidence;
import au.edu.wehi.idsv.SpanningSAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.sam.SAMRecordUtil;

import com.google.common.collect.AbstractIterator;

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
	private final Iterator<T> it;
	private boolean closed = false;
	private HashMap<String, BreakpointSummary> unpaired = new HashMap<String, BreakpointSummary>();

	public PairedEvidenceTracker(Iterator<T> it) {
		this.it = it;
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
			String partnerId = getPartnerEvidenceID((DirectedBreakpoint)evidence, evidenceId);
			if (unpaired.containsKey(evidenceId)) {
				String msg = String.format("Encountered %s multiple times.", evidenceId);
				log.error(msg);
				return false;
			}
			if (unpaired.containsKey(partnerId)) {
				BreakpointSummary remoteBp = unpaired.remove(partnerId);
				// Breakpoints must be the same
				if (!remoteBp.remoteBreakpoint().equals(((DirectedBreakpoint)evidence).getBreakendSummary())) {
					String msg = String.format("Breakpoints %s and %s differ for evidence pair %s %s", evidence.getBreakendSummary(), remoteBp, evidence.getEvidenceID(), partnerId); 
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
			String msg = String.format("Missing %d evidence pairings: ", unpaired.size()); 
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

	private static String getPartnerEvidenceID(DirectedBreakpoint evidence, String evidenceId) {
		if (evidence instanceof SpanningSAMRecordAssemblyEvidence) {
			return ((SpanningSAMRecordAssemblyEvidence)evidence).asRemote().getEvidenceID();
		} else if (evidence instanceof RealignedRemoteSoftClipEvidence || evidence instanceof RealignedRemoteSAMRecordAssemblyEvidence) {
			// remove R prefix
			if (evidenceId.charAt(0) != 'R') throw new RuntimeException("Sanity check failure: unexpected remote evidence ID " + evidenceId);
			return evidenceId.substring(1);
		} else if (evidence instanceof RealignedSoftClipEvidence || evidence instanceof RealignedSAMRecordAssemblyEvidence) {
			// add R prefix
			return "R" + evidenceId;
		} else if (evidence instanceof DiscordantReadPair) {
			String partner = swapSuffix(evidenceId, SAMRecordUtil.FIRST_OF_PAIR_NAME_SUFFIX, SAMRecordUtil.SECOND_OF_PAIR_NAME_SUFFIX);
			if (partner == null) {
				throw new RuntimeException("Sanity check failure: unexpected read pair evidence ID " + evidenceId);
			}
			return partner;
		} else if (evidence instanceof VariantContextDirectedBreakpoint) {
			String partner = swapSuffix(evidenceId, MaximalEvidenceCliqueIterator.MATE_BREAKEND_ID_SUFFIX_LOW, MaximalEvidenceCliqueIterator.MATE_BREAKEND_ID_SUFFIX_HIGH);
			if (partner == null) {
				throw new RuntimeException("Sanity check failure: unexpected VCF breakend evidence ID " + evidenceId);
			}
			return partner;
		}
		throw new RuntimeException("Sanity check failure: unhandled directed breakpoint type " + evidence.getClass().getName());
	}
	private static String swapSuffix(String evidenceId, String suffixA, String suffixB) {
		if (evidenceId.endsWith(suffixA)) {
			return evidenceId.substring(0, evidenceId.length() - suffixA.length()) + suffixB;
		} else if (evidenceId.endsWith(suffixB)) {
			return evidenceId.substring(0, evidenceId.length() - suffixB.length()) + suffixA;
		}
		return null;
	}
}
