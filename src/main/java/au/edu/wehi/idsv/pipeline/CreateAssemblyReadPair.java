package au.edu.wehi.idsv.pipeline;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ExecutorService;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.ProcessStep;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.sam.SAMFileUtil;
import au.edu.wehi.idsv.sam.SAMRecordMateCoordinateComparator;
import au.edu.wehi.idsv.util.FileHelper;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

/**
 * Creates read pairs for assembly and realignment
 * @author Daniel Cameron
 *
 */
public class CreateAssemblyReadPair extends DataTransformStep {
	private static final Log log = Log.getInstance(CreateAssemblyReadPair.class);
	private static final String UNSORTED_FILE_PREFIX = "unsorted.";
	private final AssemblyEvidenceSource source;
	private final List<SAMFileWriter> sortedWriters = Lists.newArrayList();
	private final List<SAMFileWriter> mateWriters = Lists.newArrayList();
	private final SAMFileHeader header;
	public CreateAssemblyReadPair(final ProcessingContext processContext, final AssemblyEvidenceSource source) {
		super(processContext);
		this.source = source;
		this.header = new SAMFileHeader() {{
			setSequenceDictionary(processContext.getReference().getSequenceDictionary());
		}};
	}
	@Override
	public void process(EnumSet<ProcessStep> steps) {
		process(steps, null);
	}
	public void process(EnumSet<ProcessStep> steps, ExecutorService threadpool) {
		if (isComplete() || !steps.contains(ProcessStep.SORT_REALIGNED_ASSEMBLIES)) {
			log.debug("no work to do");
			return;
		}
		if (!canProcess()) {
			String msg = String.format("Assembly realignment not completed. Unable to process");
			log.error(msg);
			throw new IllegalStateException(msg);
		}
		try {
			log.info("START: assembly remote breakend sorting ");
			createUnsortedOutputWriters();
			writeUnsortedOutput();
			closeUnsortedWriters();
			sort(threadpool);
			if (processContext.shouldProcessPerChromosome()) {
				// write out single file for debugging purposes
				assemblyPerChrToAssembly();
			}
			deleteTemp();
			log.info("SUCCESS: assembly remote breakend sorting ");
		} catch (Exception e) {
			String msg = "Unable to sort assembly breakpoints";
			log.error(e, msg);
			close();
			deleteTemp();
			deleteOutput();
			throw new RuntimeException(msg, e);
		}
	}
	private void assemblyPerChrToAssembly() {
		FileSystemContext fsc = processContext.getFileSystemContext();
		File singleFile = fsc.getAssembly(source.getFileIntermediateDirectoryBasedOn());
		SAMFileHeader header = processContext.getBasicSamHeader();
		header.setSortOrder(SortOrder.coordinate);
		SAMFileWriter writer = null;
		log.debug("Writing assemblies to single file " + singleFile);
		try {
			writer = processContext.getSamFileWriterFactory(true).makeSAMOrBAMWriter(header, true, FileSystemContext.getWorkingFileFor(singleFile));
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				File perChr = fsc.getAssemblyForChr(source.getFileIntermediateDirectoryBasedOn(), seq.getSequenceName());
				SamReader reader = null;
				try {
					reader = processContext.getSamReader(perChr);
					for (SAMRecord r : reader) {
						writer.addAlignment(r);
					}
				} finally {
					CloserUtil.close(reader);
					reader = null;
				}
			}
			writer.close();
			writer = null;
			FileHelper.move(FileSystemContext.getWorkingFileFor(singleFile), singleFile, true);
		} catch (IOException e) {
			log.debug("Error writing " + singleFile, e);
		} finally {
			CloserUtil.close(writer);
		}
	}
	private void sort(ExecutorService threadpool) throws IOException {
		FileSystemContext fsc = processContext.getFileSystemContext();
		List<SAMFileUtil.SortCallable> tasks = Lists.newArrayList();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				tasks.add(new SAMFileUtil.SortCallable(processContext,
						FileSystemContext.getWorkingFileFor((fsc.getAssemblyForChr(source.getFileIntermediateDirectoryBasedOn(), seq.getSequenceName())), UNSORTED_FILE_PREFIX),
						fsc.getAssemblyForChr(source.getFileIntermediateDirectoryBasedOn(), seq.getSequenceName()),
						SortOrder.coordinate));
				tasks.add(new SAMFileUtil.SortCallable(processContext,
						FileSystemContext.getWorkingFileFor((fsc.getAssemblyMateForChr(source.getFileIntermediateDirectoryBasedOn(), seq.getSequenceName())), UNSORTED_FILE_PREFIX),
						fsc.getAssemblyMateForChr(source.getFileIntermediateDirectoryBasedOn(), seq.getSequenceName()),
						new SAMRecordMateCoordinateComparator()));
			}
		} else {
			tasks.add(new SAMFileUtil.SortCallable(processContext,
					FileSystemContext.getWorkingFileFor((fsc.getAssembly(source.getFileIntermediateDirectoryBasedOn())), UNSORTED_FILE_PREFIX),
					fsc.getAssembly(source.getFileIntermediateDirectoryBasedOn()),
					SortOrder.coordinate));
			tasks.add(new SAMFileUtil.SortCallable(processContext,
					FileSystemContext.getWorkingFileFor((fsc.getAssemblyMate(source.getFileIntermediateDirectoryBasedOn())), UNSORTED_FILE_PREFIX),
					fsc.getAssemblyMate(source.getFileIntermediateDirectoryBasedOn()),
					new SAMRecordMateCoordinateComparator()));
		}
		//if (threadpool != null) {
		//	try {
		//		log.debug("Issuing parallel sort tasks");
		//		for (Future<Void> future : threadpool.invokeAll(tasks)) {
		//			// throw any exception
		//			future.get();
		//		}
		//		log.debug("Parallel sort tasks complete");
		//	} catch (InterruptedException e) {
		//		log.error(e);
		//		throw new RuntimeException("Interrupted waiting for VCF sort", e);
		//	} catch (ExecutionException e) {
		//		throw new RuntimeException("Failed to sort VCF", e.getCause());
		//	}
		//} else {
			log.debug("Serial sort");
			for (SAMFileUtil.SortCallable c : tasks) {
				c.call();
			}
		//}
	}
	@Override
	public void close() {
		closeUnsortedWriters();
		super.close();
	}
	private void closeUnsortedWriters() {
		for (SAMFileWriter w : Iterables.concat(sortedWriters, mateWriters)) {
			w.close();
		}
		sortedWriters.clear();
		mateWriters.clear();
	}
	private void writeUnsortedOutput() {
		Iterator<SAMRecordAssemblyEvidence> it = source.iterator(false, processContext.getAssemblyParameters().writeFilteredAssemblies);
		while (it.hasNext()) {
			SAMRecordAssemblyEvidence e = it.next();
			SAMRecord assembly = e.getSAMRecord();
			SAMRecord realign = e.getRemoteSAMRecord();
			sortedWriters.get(assembly.getReferenceIndex() % sortedWriters.size()).addAlignment(assembly);
			sortedWriters.get(realign.getReferenceIndex() % sortedWriters.size()).addAlignment(realign);
			mateWriters.get(assembly.getMateReferenceIndex() % mateWriters.size()).addAlignment(assembly);
			mateWriters.get(realign.getMateReferenceIndex() % mateWriters.size()).addAlignment(realign);
		}
	}
	private void createUnsortedOutputWriters() {
		FileSystemContext fsc = processContext.getFileSystemContext();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				sortedWriters.add(processContext.getSamFileWriterFactory(true).makeSAMOrBAMWriter(header, true,
						FileSystemContext.getWorkingFileFor((fsc.getAssemblyForChr(source.getFileIntermediateDirectoryBasedOn(), seq.getSequenceName())), UNSORTED_FILE_PREFIX)));
				mateWriters.add(processContext.getSamFileWriterFactory(false).makeSAMOrBAMWriter(header, true,
						FileSystemContext.getWorkingFileFor((fsc.getAssemblyMateForChr(source.getFileIntermediateDirectoryBasedOn(), seq.getSequenceName())), UNSORTED_FILE_PREFIX)));
			}
		} else {
			sortedWriters.add(processContext.getSamFileWriterFactory(true).makeSAMOrBAMWriter(header, true,
					FileSystemContext.getWorkingFileFor((fsc.getAssembly(source.getFileIntermediateDirectoryBasedOn())), UNSORTED_FILE_PREFIX)));
			mateWriters.add(processContext.getSamFileWriterFactory(false).makeSAMOrBAMWriter(header, true,
					FileSystemContext.getWorkingFileFor((fsc.getAssemblyMate(source.getFileIntermediateDirectoryBasedOn())), UNSORTED_FILE_PREFIX)));
		}
	}
	@Override
	protected Log getLog() {
		return log;
	}
	@Override
	public List<File> getInputs() {
		return ImmutableList.of();
	}
	@Override
	public boolean canProcess() {
		return source.isRealignmentComplete();
	}
	@Override
	public List<File> getOutput() {
		List<File> outputs = Lists.newArrayList();
		FileSystemContext fsc = processContext.getFileSystemContext();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				outputs.add(fsc.getAssemblyForChr(source.getFileIntermediateDirectoryBasedOn(), seq.getSequenceName()));
				outputs.add(fsc.getAssemblyMateForChr(source.getFileIntermediateDirectoryBasedOn(), seq.getSequenceName()));
			}
		} else {
			outputs.add(fsc.getAssembly(source.getFileIntermediateDirectoryBasedOn()));
			outputs.add(fsc.getAssemblyMate(source.getFileIntermediateDirectoryBasedOn()));
		}
		return outputs;
	}
	@Override
	public List<File> getTemporary() {
		List<File> files = Lists.newArrayList();
		FileSystemContext fsc = processContext.getFileSystemContext();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				files.add(FileSystemContext.getWorkingFileFor((fsc.getAssemblyForChr(source.getFileIntermediateDirectoryBasedOn(), seq.getSequenceName())), UNSORTED_FILE_PREFIX));
				files.add(FileSystemContext.getWorkingFileFor((fsc.getAssemblyMateForChr(source.getFileIntermediateDirectoryBasedOn(), seq.getSequenceName())), UNSORTED_FILE_PREFIX));
				files.add(FileSystemContext.getWorkingFileFor((fsc.getAssemblyForChr(source.getFileIntermediateDirectoryBasedOn(), seq.getSequenceName()))));
				files.add(FileSystemContext.getWorkingFileFor((fsc.getAssemblyMateForChr(source.getFileIntermediateDirectoryBasedOn(), seq.getSequenceName()))));
			}
		} else {
			files.add(FileSystemContext.getWorkingFileFor((fsc.getAssembly(source.getFileIntermediateDirectoryBasedOn())), UNSORTED_FILE_PREFIX));
			files.add(FileSystemContext.getWorkingFileFor((fsc.getAssemblyMate(source.getFileIntermediateDirectoryBasedOn())), UNSORTED_FILE_PREFIX));
			files.add(FileSystemContext.getWorkingFileFor((fsc.getAssembly(source.getFileIntermediateDirectoryBasedOn()))));
			files.add(FileSystemContext.getWorkingFileFor((fsc.getAssemblyMate(source.getFileIntermediateDirectoryBasedOn()))));
		}
		return files;
	}
}
