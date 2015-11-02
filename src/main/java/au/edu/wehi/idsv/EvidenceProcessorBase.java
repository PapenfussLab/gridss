package au.edu.wehi.idsv;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;

import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.primitives.Ints;

import au.edu.wehi.idsv.util.AutoClosingMergedIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;

public abstract class EvidenceProcessorBase {
	private static final Log log = Log.getInstance(EvidenceProcessorBase.class);
	protected final ProcessingContext processContext;
	protected final File output;
	protected final List<SAMEvidenceSource> samEvidence;
	protected final AssemblyEvidenceSource assemblyEvidence;
	protected final List<Closeable> toClose = Lists.newArrayList();
	public EvidenceProcessorBase(ProcessingContext context, File output, List<SAMEvidenceSource> samEvidence, AssemblyEvidenceSource assemblyEvidence) {
		this.processContext = context;
		this.output = output;
		this.samEvidence = samEvidence;
		this.assemblyEvidence = assemblyEvidence;
	}
	@SuppressWarnings("unchecked")
	protected CloseableIterator<DirectedEvidence> getAllEvidence(
			boolean includeAssembly,
			boolean includeAssemblyRemote,
			boolean includeReadPair,
			boolean includeSoftClip,
			boolean includeSoftClipRemote,
			boolean includeFilteredAssemblies) {
		List<CloseableIterator<DirectedEvidence>> evidenceItList = Lists.newArrayList();
		if (includeAssembly && assemblyEvidence != null) {
			evidenceItList.add((CloseableIterator<DirectedEvidence>)(Object)assemblyEvidence.iterator(includeAssemblyRemote, includeFilteredAssemblies));
		}
		if (includeReadPair || includeSoftClip || includeSoftClipRemote) {
			for (SAMEvidenceSource source : samEvidence) {
				evidenceItList.add(source.iterator(includeReadPair, includeSoftClip, includeSoftClipRemote));
			}
		}
		toClose.addAll(evidenceItList);
		CloseableIterator<DirectedEvidence> evidenceIt = new AutoClosingMergedIterator<DirectedEvidence>(evidenceItList, DirectedEvidenceOrder.ByNatural);
		return evidenceIt;
	}
	protected CloseableIterator<DirectedEvidence> getAllEvidence(
			boolean includeAssembly,
			boolean includeAssemblyRemote,
			boolean includeReadPair,
			boolean includeSoftClip,
			boolean includeSoftClipRemote) {
		return getAllEvidence(includeAssembly, includeAssemblyRemote, includeReadPair, includeSoftClip, includeSoftClipRemote, false);
	}
	@SuppressWarnings("unchecked")
	protected CloseableIterator<DirectedEvidence> getEvidenceForChr(
			boolean includeAssembly,
			boolean includeAssemblyRemote,
			boolean includeReadPair,
			boolean includeSoftClip,
			boolean includeSoftClipRemote,
			int... referenceIndex) {
		SortedSet<Integer> refList = Sets.newTreeSet(Ints.asList(referenceIndex));
		List<Iterator<DirectedEvidence>> evidenceItList = Lists.newArrayList();
		if (includeAssembly && assemblyEvidence != null) {
			List<CloseableIterator<SAMRecordAssemblyEvidence>> itList = Lists.newArrayList();
			for (int i : refList) {
				String chr = processContext.getReference().getSequenceDictionary().getSequence(i).getSequenceName();
				CloseableIterator<SAMRecordAssemblyEvidence> assIt = assemblyEvidence.iterator(includeAssemblyRemote, false, chr);
				toClose.add(assIt);
				itList.add(assIt);
			}
			evidenceItList.add((Iterator<DirectedEvidence>)(Object)Iterators.concat(itList.iterator()));
		}
		if (includeReadPair || includeSoftClip || includeSoftClipRemote) {
			for (SAMEvidenceSource source : samEvidence) {
				for (int i : refList) {
					String chr = processContext.getReference().getSequenceDictionary().getSequence(i).getSequenceName();
					CloseableIterator<DirectedEvidence> samIt = source.iterator(includeReadPair, includeSoftClip, includeSoftClipRemote, chr);
					evidenceItList.add(samIt);
				}
			}
		}
		CloseableIterator<DirectedEvidence> evidenceIt = new AutoClosingMergedIterator<DirectedEvidence>(evidenceItList, DirectedEvidenceOrder.ByNatural);
		return evidenceIt;
	}
	public void close() {
		for (Closeable c : toClose) {
			if (c != null) {
				try {
					c.close();
				} catch (IOException e) {
					log.error(e, " error closing ", c);
				}
			}
		}
		toClose.clear();
	}
	public abstract void process();
}