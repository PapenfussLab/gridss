package au.edu.wehi.idsv;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

/**
 * Maps reads to their associated assemblies.
 * 
 * Note: this iterator requires that the exact evidence supplied to the
 * original assemblies is supplied to this iterator. Failure to do so will
 * result in the missing evidenceIDs being retained in memory.  
 * 
 * @author Daniel Cameron
 *
 */
public class AssemblyAssociator implements CloseableIterator<DirectedEvidence> {
	private final Iterator<DirectedEvidence> it;
	private final Iterator<SAMRecord> assit;
	private final int windowSize;
	private final HashMap<String, String> evidenceToAssemblyName = new HashMap<>();
	private SAMRecord lastAssembly = null;;
	public AssemblyAssociator(Iterator<DirectedEvidence> it, Iterator<SAMRecord> rawAssemblies, int windowSize) {
		this.it = it;
		this.assit = rawAssemblies;
		this.windowSize = windowSize;
	}
	@Override
	public boolean hasNext() {
		return it.hasNext();
	}

	@Override
	public DirectedEvidence next() {
		DirectedEvidence e = it.next();
		associate(e);
		return e;
	}
	private DirectedEvidence associate(DirectedEvidence e) {
		ensureAssembliesLoadedUntil(e.getBreakendSummary());
		setAssociatedAssembly(e, evidenceToAssemblyName.remove(e.getEvidenceID()));
		flushBefore(e.getBreakendSummary());
		return e;
	}
	private void flushBefore(BreakendSummary breakendSummary) {
		// See class description. No flushing is required if all evidence is passed in.
	}
	private void ensureAssembliesLoadedUntil(BreakendSummary breakendSummary) {
		while (assit.hasNext() && (lastAssembly == null || !isAfter(breakendSummary, lastAssembly))) {
			lastAssembly = assit.next();
			load(lastAssembly);
		}
	}
	private void load(SAMRecord ass) {
		assert(ass != null);
		Collection<String> eids = new AssemblyAttributes(ass).getEvidenceIDs();
		for (String eid : eids) {
			evidenceToAssemblyName.put(eid, ass.getReadName());
		}
	}
	private boolean isAfter(BreakendSummary breakendSummary, SAMRecord position) {
		return position.getReferenceIndex() > breakendSummary.referenceIndex ||
			(position.getReferenceIndex() == breakendSummary.referenceIndex && position.getUnclippedStart() > breakendSummary.end + windowSize);
	}
	@Override
	public void close() {
		CloserUtil.close(it);
		CloserUtil.close(assit);
	}
	private void setAssociatedAssembly(DirectedEvidence e, String assemblyName) {
		if (e instanceof SingleReadEvidence) {
			((SingleReadEvidence)e).setAssociatedAssemblyName(assemblyName);
		} else if (e instanceof NonReferenceReadPair) {
			((NonReferenceReadPair)e).setAssociatedAssemblyName(assemblyName);
		}
	}
}
