package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import au.edu.wehi.idsv.configuration.RealignmentConfiguration;
import au.edu.wehi.idsv.visualisation.TrackedBuffer;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;

public class SAMRecordAssemblyEvidenceIterator extends AbstractIterator<SAMRecordAssemblyEvidence> implements CloseableIterator<SAMRecordAssemblyEvidence>, TrackedBuffer {
	private final ProcessingContext processContext;
	private final AssemblyEvidenceSource source;
	private final Iterator<SAMRecord> it;
	private final List<? extends Iterator<SAMRecord>> rit;
	private final List<SequentialRealignedBreakpointFactory> factories;
	private boolean includeBothBreakendsOfSpanningAssemblies;
	public SAMRecordAssemblyEvidenceIterator(
			ProcessingContext processContext,
			AssemblyEvidenceSource source,
			Iterator<SAMRecord> it,
			List<? extends Iterator<SAMRecord>> realignedIt,
			boolean includeBothBreakendsOfSpanningAssemblies) {
		this.processContext = processContext;
		this.source = source;
		this.it = it;
		this.rit = realignedIt;
		this.factories = realignedIt.stream().map(x -> new SequentialRealignedBreakpointFactory(Iterators.peekingIterator(x))).collect(Collectors.toList());
		this.includeBothBreakendsOfSpanningAssemblies = includeBothBreakendsOfSpanningAssemblies;
	}
	private SAMRecordAssemblyEvidence buffer = null;
	@Override
	protected SAMRecordAssemblyEvidence computeNext() {
		if (buffer != null) {
			SAMRecordAssemblyEvidence r = buffer;
			buffer = null;
			return r;
		}
		while (it.hasNext()) {
			SAMRecord record = it.next();
			SAMRecordAssemblyEvidence evidence = AssemblyFactory.hydrate(source, record);
			if (!(evidence instanceof DirectedBreakpoint)) {
				RealignmentConfiguration rp = evidence.getEvidenceSource().getContext().getRealignmentParameters();
				List<SAMRecord> realignments = new ArrayList<SAMRecord>(source.getRealignmentIterationCount());
				for (SequentialRealignedBreakpointFactory factory : factories) {
					Set<SAMRecord> realigned = factory.findAllAssociatedSAMRecords(evidence,
							rp.requireRealignment && 
							rp.shouldRealignBreakend(evidence) &&
							factory == factories.get(0)); // only primary realignment is required
					if (realigned != null) {
						realignments.addAll(realigned);
					}
				}
				if (realignments.size() > 0) {
					evidence = AssemblyFactory.incorporateRealignment(processContext, evidence, realignments);
				}
			}
			if (includeBothBreakendsOfSpanningAssemblies && evidence.isSpanningAssembly()) {
				buffer = ((SmallIndelSAMRecordAssemblyEvidence)evidence).asRemote();
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
		for (Iterator<SAMRecord> i : rit) {
			CloserUtil.close(i);
		}
	}
	@Override
	public void setTrackedBufferContext(String context) {
		factories.stream().forEach(factory -> factory.setTrackedBufferContext(context));
	}
	@Override
	public List<NamedTrackedBuffer> currentTrackedBufferSizes() {
		return factories.stream().flatMap(factory -> factory.currentTrackedBufferSizes().stream()).collect(Collectors.toList());
	}
}
