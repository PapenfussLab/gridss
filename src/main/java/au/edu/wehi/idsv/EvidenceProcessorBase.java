package au.edu.wehi.idsv;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;

import au.edu.wehi.idsv.util.AutoClosingMergedIterator;

import com.google.common.collect.Iterators;
import com.google.common.collect.Sets;
import com.google.common.primitives.Ints;

public abstract class EvidenceProcessorBase implements Closeable {
	private static final Log log = Log.getInstance(EvidenceProcessorBase.class);
	protected final ProcessingContext processContext;
	protected final File output;
	protected final List<EvidenceSource> evidence;
	protected final List<Closeable> toClose = new ArrayList<>();
	public EvidenceProcessorBase(ProcessingContext context, File output, List<EvidenceSource> evidence) {
		this.processContext = context;
		this.output = output;
		this.evidence = evidence;
	}
	protected CloseableIterator<DirectedEvidence> getAllEvidence(boolean assemblyOnly, boolean includeRemoteAssembly) {
		List<Iterator<DirectedEvidence>> evidenceItList = new ArrayList<>();
		for (EvidenceSource source : evidence) {
			if (source instanceof AssemblyEvidenceSource) {
				((AssemblyEvidenceSource)source).setIncludeRemote(includeRemoteAssembly);
			} else if (assemblyOnly) {
				// only include assembly evidence if assemblyOnly is set
				continue;
			}
			Iterator<DirectedEvidence> it = source.iterator();
			if (it instanceof Closeable) toClose.add((Closeable)it);
			evidenceItList.add(it);
		}
		CloseableIterator<DirectedEvidence> evidenceIt = new AutoClosingMergedIterator<DirectedEvidence>(evidenceItList, DirectedEvidenceOrder.ByNatural);
		return evidenceIt;
	}
	protected CloseableIterator<DirectedEvidence> getEvidenceForChr(boolean assemblyOnly, int... referenceIndex) {
		SortedSet<Integer> refList = Sets.newTreeSet(Ints.asList(referenceIndex));
		List<Iterator<DirectedEvidence>> evidenceItList = new ArrayList<>();
		for (EvidenceSource source : evidence) {
			if (assemblyOnly && !(source instanceof AssemblyEvidenceSource)) continue;
			List<Iterator<DirectedEvidence>> sourceIt = new ArrayList<>();
			for (int i : refList) {
				CloseableIterator<DirectedEvidence> it = source.iterator(processContext.getReference().getSequenceDictionary().getSequence(i).getSequenceName());
				toClose.add(it);
				sourceIt.add(it);
			}
			evidenceItList.add(Iterators.concat(sourceIt.iterator()));
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