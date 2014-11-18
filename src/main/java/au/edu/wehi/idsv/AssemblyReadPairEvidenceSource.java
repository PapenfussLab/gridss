package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;

import java.io.File;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ExecutorService;

import au.edu.wehi.idsv.pipeline.CreateAssemblyReadPair;

import com.google.common.base.Function;

/**
 * Assembly evidence source that iterates over both local and remote assembly evidence
 * using proxy SAMRecords for evidence
 * @author Daniel Cameron
 *
 */
public class AssemblyReadPairEvidenceSource extends AssemblyEvidenceSource {
	public AssemblyReadPairEvidenceSource(ProcessingContext processContext, List<SAMEvidenceSource> evidence, File intermediateFileLocation) {
		super(processContext, evidence, intermediateFileLocation);
	}
	@Override
	public void ensureAssembled(ExecutorService threadpool) {
		super.ensureAssembled(threadpool);
		if (isRealignmentComplete()) {
			CreateAssemblyReadPair step = new CreateAssemblyReadPair(processContext, this);
			step.process(EnumSet.of(ProcessStep.SORT_REALIGNED_ASSEMBLIES), threadpool);
			step.close();
		}
	}
	public CloseableIterator<SAMRecordAssemblyEvidence> iterator(final boolean includeFiltered) {
		return iterator(false, includeFiltered);
	}
	public CloseableIterator<SAMRecordAssemblyEvidence> iterator(final boolean includeRemote, final boolean includeFiltered) {
		if (processContext.shouldProcessPerChromosome()) {
			// Lazily iterator over each input
			return new PerChromosomeAggregateIterator<SAMRecordAssemblyEvidence>(processContext.getReference().getSequenceDictionary(), new Function<String, Iterator<SAMRecordAssemblyEvidence>>() {
				@Override
				public Iterator<SAMRecordAssemblyEvidence> apply(String chr) {
					return perChrIterator(includeRemote, includeFiltered, chr);
				}
			});
		} else {
			return singleFileIterator(includeRemote, includeFiltered);
		}
	}
	public CloseableIterator<SAMRecordAssemblyEvidence> iterator(boolean includeRemote, boolean includeFiltered, String chr) {
		if (processContext.shouldProcessPerChromosome()) {
			return perChrIterator(includeRemote, includeFiltered, chr);
		} else {
			return new ChromosomeFilteringIterator<SAMRecordAssemblyEvidence>(singleFileIterator(includeRemote, includeFiltered), processContext.getDictionary().getSequence(chr).getSequenceIndex(), true);
		}
	}
	protected CloseableIterator<SAMRecordAssemblyEvidence> perChrIterator(boolean includeRemote, boolean includeFiltered, String chr) {
		FileSystemContext fsc = processContext.getFileSystemContext();
		return samAssemblyRealignIterator(
				includeRemote,
				includeFiltered,
				fsc.getAssemblyForChr(input, chr),
				fsc.getAssemblyMateForChr(input, chr),
				chr);
	}
	protected CloseableIterator<SAMRecordAssemblyEvidence> singleFileIterator(boolean includeRemote, boolean includeFiltered) {
		FileSystemContext fsc = processContext.getFileSystemContext();
		return samAssemblyRealignIterator(
				includeRemote,
				includeFiltered,
				fsc.getAssembly(input),
				fsc.getAssemblyMate(input),
				"");
	}
	private CloseableIterator<SAMRecordAssemblyEvidence> samAssemblyRealignIterator(
			boolean includeRemote,
			boolean includeFiltered,
			File sorted,
			File mateSorted,
			String chr) {
		CloseableIterator<SAMRecord> it = processContext.getSamReaderIterator(processContext.getFileSystemContext().getAssembly(getFileIntermediateDirectoryBasedOn()));
		CloseableIterator<SAMRecord> mateIt = processContext.getSamReaderIterator(processContext.getFileSystemContext().getAssemblyMate(getFileIntermediateDirectoryBasedOn()));
		return new SAMRecordAssemblyEvidenceReadPairIterator(processContext, this, it, mateIt, includeRemote, includeFiltered);
	}
}
