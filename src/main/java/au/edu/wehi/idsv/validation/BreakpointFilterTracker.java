package au.edu.wehi.idsv.validation;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import java.util.Collection;
import java.util.Iterator;
import java.util.Set;

import au.edu.wehi.idsv.DirectedBreakpoint;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;

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
	private final Multimap<String, String> filters = TreeMultimap.create();
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
		if (!it.hasNext()) return endOfData();
		T evidence = it.next();
		assert(trackFilters(evidence));
		return evidence;
	}

	private boolean trackFilters(T evidence) {
		if (evidence instanceof DirectedBreakpoint) {
			String eventid = evidence.getAttributeAsString(VcfSvConstants.BREAKEND_EVENT_ID_KEY, "");
			Set<String> currentFilters = evidence.getFilters();
			if (filters.containsKey(eventid)) {
				Collection<String> storedFiltered = filters.get(eventid);
				if (currentFilters.size() != storedFiltered.size() || !currentFilters.containsAll(storedFiltered)) {
					String msg = "Breakend filter mismatch for " + eventid + " " + toString(storedFiltered) + " vs " + toString(currentFilters); 
					log.error(msg);
					return !doAssert;
				}
				filters.removeAll(eventid);
			} else {
				filters.putAll(eventid, currentFilters);
			}
		}
		return true;
	}
	private static String toString(Collection<String> filters) {
		StringBuilder sb = new StringBuilder("{");
		for (String s : filters) {
			sb.append(s);
			sb.append(", ");
		}
		sb.append("}");
		return sb.toString();
	}
}
