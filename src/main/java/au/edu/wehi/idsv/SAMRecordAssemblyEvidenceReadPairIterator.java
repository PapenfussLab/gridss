package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;

public class SAMRecordAssemblyEvidenceReadPairIterator extends AbstractIterator<SAMRecordAssemblyEvidence> implements CloseableIterator<SAMRecordAssemblyEvidence> {
	private final ProcessingContext processContext;
	private final AssemblyEvidenceSource source;
	private final CloseableIterator<SAMRecord> it;
	private final CloseableIterator<SAMRecord> rit;
	private final SequentialAssemblyReadPairFactory factory;
	private final boolean includeRemote;
	private final boolean includeFiltered;
	public SAMRecordAssemblyEvidenceReadPairIterator(
			ProcessingContext processContext,
			AssemblyEvidenceSource source,
			CloseableIterator<SAMRecord> it,
			CloseableIterator<SAMRecord> mateIt,
			boolean includeRemote,
			boolean includeFiltered) {
		this.processContext = processContext;
		this.source = source;
		this.it = it;
		this.rit = mateIt;
		this.factory = new SequentialAssemblyReadPairFactory(Iterators.peekingIterator(mateIt));
		this.includeFiltered = includeFiltered;
		this.includeRemote = includeRemote;
	}
	@Override
	protected SAMRecordAssemblyEvidence computeNext() {
		while (it.hasNext()) {
			SAMRecord record = it.next();
			if (record.getFirstOfPairFlag() || includeRemote) {
				SAMRecordAssemblyEvidence evidence = factory.createAssembly(record, processContext, source);
				if (evidence != null) {
					processContext.getAssemblyParameters().applyFilters(evidence);
					if (!evidence.isAssemblyFiltered() || includeFiltered) {
						return evidence;
					}
				}
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