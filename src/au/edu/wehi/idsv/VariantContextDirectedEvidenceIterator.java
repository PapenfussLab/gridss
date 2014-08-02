package au.edu.wehi.idsv;

import java.util.Iterator;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;

public class VariantContextDirectedEvidenceIterator extends AbstractIterator<VariantContextDirectedEvidence> implements CloseableIterator<VariantContextDirectedEvidence> {
	private final ProcessingContext processContext;
	private final EvidenceSource source;
	private final Iterator<VariantContext> it;
	private final Iterator<SAMRecord> rit;
	private final SequentialRealignedBreakpointFactory factory;
	public VariantContextDirectedEvidenceIterator(
			ProcessingContext processContext,
			EvidenceSource source,
			Iterator<VariantContext> it,
			Iterator<SAMRecord> realignedIt) {
		this.processContext = processContext;
		this.source = source;
		this.it = it;
		this.rit = realignedIt;
		this.factory = new SequentialRealignedBreakpointFactory(Iterators.peekingIterator(realignedIt));
	}
	@Override
	protected VariantContextDirectedEvidence computeNext() {
		while (it.hasNext()) {
			VariantContext variant = it.next();
			IdsvVariantContext managedContext = IdsvVariantContext.create(processContext, source, variant);
			if (managedContext instanceof DirectedEvidence) {
				VariantContextDirectedEvidence assembly = (VariantContextDirectedEvidence)variant;
				SAMRecord realigned = factory.findAssociatedSAMRecord(assembly);
				assembly = AssemblyBuilder.incorporateRealignment(processContext, assembly, realigned);
				return assembly;
			}
		}
		return endOfData();
	}
	@Override
	public void close() {
		CloserUtil.close(it);
		CloserUtil.close(rit);
	}
}
