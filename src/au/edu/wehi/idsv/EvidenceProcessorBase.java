package au.edu.wehi.idsv;

import htsjdk.samtools.util.Log;

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

public abstract class EvidenceProcessorBase implements Closeable {
	private static final Log log = Log.getInstance(EvidenceProcessorBase.class);
	protected final ProcessingContext processContext;
	protected final File output;
	protected final List<EvidenceSource> evidence;
	protected final List<Closeable> toClose = Lists.newArrayList();
	public EvidenceProcessorBase(ProcessingContext context, File output, List<EvidenceSource> evidence) {
		this.processContext = context;
		this.output = output;
		this.evidence = evidence;
	}
	protected Iterator<DirectedEvidence> getAllEvidence() {
		List<Iterator<DirectedEvidence>> evidenceItList = Lists.newArrayList();
		for (EvidenceSource source : evidence) {
			Iterator<DirectedEvidence> it = source.iterator();
			if (it instanceof Closeable) toClose.add((Closeable)it);
			evidenceItList.add(it);
		}
		Iterator<DirectedEvidence> evidenceIt = Iterators.mergeSorted(evidenceItList, DirectedEvidenceOrder.ByNatural);
		return evidenceIt;
	}
	protected Iterator<DirectedEvidence> getEvidenceForChr(int... referenceIndex) {
		SortedSet<Integer> refList = Sets.newTreeSet(Ints.asList(referenceIndex));
		List<Iterator<DirectedEvidence>> evidenceItList = Lists.newArrayList();
		for (EvidenceSource source : evidence) {
			List<Iterator<DirectedEvidence>> sourceIt = Lists.newArrayList();
			for (int i : refList) {
				Iterator<DirectedEvidence> it = source.iterator(processContext.getReference().getSequenceDictionary().getSequence(i).getSequenceName());
				if (it instanceof Closeable) toClose.add((Closeable)it);
				sourceIt.add(it);
			}
			evidenceItList.add(Iterators.concat(sourceIt.iterator()));
		}
		Iterator<DirectedEvidence> evidenceIt = Iterators.mergeSorted(evidenceItList, DirectedEvidenceOrder.ByNatural);
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