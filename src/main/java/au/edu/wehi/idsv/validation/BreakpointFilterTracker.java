package au.edu.wehi.idsv.validation;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import java.util.Iterator;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;

import au.edu.wehi.idsv.DirectedBreakpoint;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Maps;

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
		if (be1.getFilters().equals(be2.getFilters()) &&
				breakendsMatchs(be1, be2) &&
				Math.abs(be1.getPhredScaledQual() - be2.getPhredScaledQual()) <= 0.01
				
						) {
			return true;
		}
		String msg = String.format("Breakend mismatch between %s (%s){%s} and %s (%s){%s}",
				be1.getID(), be1.getBreakendSummary(), be1.getFilters(),
				be2.getID(), be2.getBreakendSummary(), be2.getFilters());
		log.error(msg);
		return !doAssert;
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

