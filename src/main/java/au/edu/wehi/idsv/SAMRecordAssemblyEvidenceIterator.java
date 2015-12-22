package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Queue;
import java.util.Set;
import java.util.stream.Collectors;

import au.edu.wehi.idsv.configuration.RealignmentConfiguration;
import au.edu.wehi.idsv.visualisation.TrackedBuffer;

import com.google.common.collect.Iterators;

public class SAMRecordAssemblyEvidenceIterator implements
		CloseableIterator<SAMRecordAssemblyEvidence>, TrackedBuffer {
	private final AssemblyEvidenceSource source;
	private final Iterator<SAMRecord> it;
	private final List<? extends Iterator<SAMRecord>> rit;
	private final List<SequentialRealignedBreakpointFactory> factories;
	private final Queue<SAMRecordAssemblyEvidence> buffer = new ArrayDeque<SAMRecordAssemblyEvidence>();
	private final boolean includeSpannedIndels;
	private final boolean includeIndelParents;
	public SAMRecordAssemblyEvidenceIterator(
			AssemblyEvidenceSource source, Iterator<SAMRecord> it,
			List<? extends Iterator<SAMRecord>> realignedIt,
			boolean includeSpannedIndels, boolean includeIndelParents) {
		this.source = source;
		this.it = it;
		this.rit = realignedIt;
		this.factories = realignedIt
				.stream()
				.map(x -> new SequentialRealignedBreakpointFactory(Iterators
						.peekingIterator(x))).collect(Collectors.toList());
		this.includeSpannedIndels = includeSpannedIndels;
		this.includeIndelParents = includeIndelParents;
	}

	private void ensureBuffer() {
		while (buffer.isEmpty() && it.hasNext()) {
			process(it.next());
		}
	}
	
	private void process(SAMRecord record) {
		SAMRecordAssemblyEvidence evidence = AssemblyFactory.hydrate(source, record);
		assert(evidence != null);
		if (evidence.getBreakendSummary() != null) {
			RealignmentConfiguration rp = source.getContext().getRealignmentParameters();
			List<SAMRecord> realignments = new ArrayList<SAMRecord>(source.getRealignmentIterationCount());
			for (SequentialRealignedBreakpointFactory factory : factories) {
				// only primary realignment is required
				boolean expectAlignment = rp.shouldRealignBreakend(evidence) && factory == factories.get(0);
				Set<SAMRecord> realigned = factory.findAllAssociatedSAMRecords(evidence, expectAlignment);
				if (realigned != null) {
					realignments.addAll(realigned);
				}
			}
			if (realignments.size() > 0) {
				evidence = AssemblyFactory.incorporateRealignment(source.getContext(), evidence, realignments);
			}
		}
		List<SpanningSAMRecordAssemblyEvidence> indels = null;
		if (includeSpannedIndels) {
			indels = evidence.getSpannedIndels();
			buffer.addAll(indels);
		}
		if (!evidence.isReferenceAssembly() && (evidence.getBreakendSummary() != null || (includeIndelParents && !indels.isEmpty()))) {
			buffer.add(evidence);
		}
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
		factories.stream().forEach(
				factory -> factory.setTrackedBufferContext(context));
	}

	@Override
	public List<NamedTrackedBuffer> currentTrackedBufferSizes() {
		return factories
				.stream()
				.flatMap(
						factory -> factory.currentTrackedBufferSizes().stream())
				.collect(Collectors.toList());
	}

	@Override
	public boolean hasNext() {
		ensureBuffer();
		return !buffer.isEmpty();
	}

	@Override
	public SAMRecordAssemblyEvidence next() {
		ensureBuffer();
		return buffer.poll();
	}
}
