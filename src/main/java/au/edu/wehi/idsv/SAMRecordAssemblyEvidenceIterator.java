package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.util.Iterator;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;

public class SAMRecordAssemblyEvidenceIterator extends AbstractIterator<SAMRecordAssemblyEvidence> implements CloseableIterator<SAMRecordAssemblyEvidence> {
	private final ProcessingContext processContext;
	private final AssemblyEvidenceSource source;
	private final Iterator<SAMRecord> it;
	private final Iterator<SAMRecord> rit;
	private final SequentialRealignedBreakpointFactory factory;
	public SAMRecordAssemblyEvidenceIterator(
			ProcessingContext processContext,
			AssemblyEvidenceSource source,
			Iterator<SAMRecord> it,
			Iterator<SAMRecord> realignedIt) {
		this.processContext = processContext;
		this.source = source;
		this.it = it;
		this.rit = realignedIt;
		this.factory = realignedIt != null ? new SequentialRealignedBreakpointFactory(Iterators.peekingIterator(this.rit)) : null;
	}
	@Override
	protected SAMRecordAssemblyEvidence computeNext() {
		while (it.hasNext()) {
			SAMRecord record = it.next();
			SAMRecordAssemblyEvidence evidence = AssemblyFactory.hydrate(source, record);
			if (factory != null && !(evidence instanceof DirectedBreakpoint)) {
				RealignmentParameters rp = evidence.getEvidenceSource().getContext().getRealignmentParameters();
				SAMRecord realigned = factory.findAssociatedSAMRecord(evidence,
						rp.requireRealignment && 
						rp.shouldRealignBreakend(evidence));
				evidence = AssemblyFactory.incorporateRealignment(processContext, evidence, realigned);
			}
			if (evidence != null && !evidence.isReferenceAssembly()) {
				return evidence;
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
