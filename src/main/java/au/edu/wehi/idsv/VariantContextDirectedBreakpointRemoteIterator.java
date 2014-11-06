package au.edu.wehi.idsv;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.Iterator;

public class VariantContextDirectedBreakpointRemoteIterator implements CloseableIterator<VariantContextDirectedBreakpointRemote> {
	private final ProcessingContext processContext;
	private final EvidenceSource source;
	private final Iterator<? extends VariantContext> it;
	private boolean closed = false;
	public VariantContextDirectedBreakpointRemoteIterator(
			ProcessingContext processContext,
			EvidenceSource source,
			Iterator<? extends VariantContext> it) {
		this.processContext = processContext;
		this.source = source;
		this.it = it;
	}
	@Override
	public void close() {
		CloserUtil.close(it);
	}
	@Override
	public boolean hasNext() {
		return !closed && it.hasNext();
	}
	@Override
	public VariantContextDirectedBreakpointRemote next() {
		return new VariantContextDirectedBreakpointRemote(processContext, source, it.next());
	}
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
}
