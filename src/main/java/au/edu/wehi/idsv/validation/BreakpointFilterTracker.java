package au.edu.wehi.idsv.validation;

import java.util.Iterator;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Maps;

import au.edu.wehi.idsv.DirectedBreakpoint;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

/**
 * Pass-through iterator that tracks filters applied to breakpoints 
 * to ensure filters at both breakends match. 
 * 
 * Tracking is performed via assertion. If assertions are not enabled,
 * no tracking is performed.
 * 
 * @author Daniel Cameron
 *
 */
public class BreakpointFilterTracker<T extends VariantContextDirectedEvidence> extends AbstractIterator<T> implements CloseableIterator<T> {
	private static final Log log = Log.getInstance(BreakpointFilterTracker.class);
	private final Iterator<T> it;
	private final boolean doAssert;
	private final Map<String, T> partner = Maps.newHashMap();
	private boolean closed = false;
	public BreakpointFilterTracker(Iterator<T> it, boolean doAssert) {
		this.it = it;
		this.doAssert = doAssert;
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
		if (!it.hasNext()) {
			assert(allMatched());
			return endOfData();
		}
		T evidence = it.next();
		assert(track(evidence));
		return evidence;
	}

	private boolean track(T evidence) {
		String mateid = evidence.getAttributeAsString(VcfSvConstants.MATE_BREAKEND_ID_KEY, null);
		if (!StringUtils.isEmpty(mateid)) {
			if (partner.containsKey(mateid)) {
				return isValidBreakpoint(evidence, partner.remove(mateid));
			} else {
				partner.put(evidence.getID(), evidence);
			}
		}
		return true;
	}
	private boolean isValidBreakpoint(T be1, T be2) {
		String be1value = "";
		String be2value = "";
		if (!be1.getFilters().equals(be2.getFilters())) {
			be1value += " {" + be1.getFilters().toString() + "}";
			be2value += " {" + be2.getFilters().toString() + "}";
		}
		if (!breakendsMatchs(be1, be2)) {
			be1value += be1.getBreakendSummary().toString();
			be2value += be2.getBreakendSummary().toString();
		}
		if (Math.abs(be1.getPhredScaledQual() - be2.getPhredScaledQual()) > 0.01) {
			be1value += " Q:" + be1.getPhredScaledQual();
			be2value += " Q:" + be2.getPhredScaledQual();
		}
		if (be1 instanceof VariantContextDirectedBreakpoint && be2 instanceof VariantContextDirectedBreakpoint) {
			VariantContextDirectedBreakpoint bp1 = (VariantContextDirectedBreakpoint) be1;
			VariantContextDirectedBreakpoint bp2 = (VariantContextDirectedBreakpoint) be2;
			if (bp1.getBreakpointEvidenceCountAssembly() != bp2.getBreakpointEvidenceCountAssembly()) {
				be1value += " AS:" + bp1.getBreakpointEvidenceCountLocalAssembly() + " RAS:" + bp1.getBreakpointEvidenceCountRemoteAssembly();
				be2value += " AS:" + bp2.getBreakpointEvidenceCountLocalAssembly() + " RAS:" + bp2.getBreakpointEvidenceCountRemoteAssembly();
			}
			if (bp1.getBreakpointEvidenceCountSoftClip() != bp2.getBreakpointEvidenceCountSoftClip()) {
				be1value += " SC:" + bp1.getBreakpointEvidenceCountLocalSoftClip() + " RSC:" + bp1.getBreakpointEvidenceCountRemoteSoftClip();
				be2value += " SC:" + bp2.getBreakpointEvidenceCountLocalSoftClip() + " RSC:" + bp2.getBreakpointEvidenceCountRemoteSoftClip();
			}
			if (bp1.getBreakpointEvidenceCountReadPair() != bp2.getBreakpointEvidenceCountReadPair()) {
				be1value += " RP:" + bp1.getBreakpointEvidenceCountReadPair();
				be2value += " RP:" + bp2.getBreakpointEvidenceCountReadPair();
			}
			if (bp1.getHomologySequence().length() != bp2.getHomologySequence().length()) {
				be1value += " HOMSEQ:" + bp1.getHomologySequence();
				be2value += " HOMSEQ:" + bp2.getHomologySequence();
			}
		}
		if (be1value.length() > 0) {
			String msg = String.format("Breakend mismatch between %s %s and %s %s",
					be1.getID(), be1value,
					be2.getID(), be2value);
			log.error(msg);
			return !doAssert;
		}
		return true;
	}
	private boolean breakendsMatchs(T be1, T be2) {
		if (be1 instanceof DirectedBreakpoint && be2 instanceof DirectedBreakpoint) {
			return ((DirectedBreakpoint)be1).getBreakendSummary().remoteBreakpoint().equals(((DirectedBreakpoint)be2).getBreakendSummary());
		}
		return true;
	}
	private boolean allMatched() {
		if (!partner.isEmpty()) {
			String msg = String.format("Missing %d breakend partners: ", partner.size()); 
			StringBuilder sb = new StringBuilder();
			sb.append(msg);
			for (String s : partner.keySet()) {
				sb.append(s);
				sb.append(", ");
			}
			log.error(sb.toString());
			return false;
		}
		return true;
	}
}

