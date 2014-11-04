package au.edu.wehi.idsv.pipeline;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import java.io.File;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.ProcessStep;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.vcf.VcfFileUtil.SortCallable;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;

/**
 * Sorts assembly breakends by the remote breakend coordinate
 * @author Daniel Cameron
 *
 */
public class SortRealignedAssemblies extends DataTransformStep {
	private static final Log log = Log.getInstance(SortRealignedAssemblies.class);
	private final AssemblyEvidenceSource source;
	private final List<VariantContextWriter> asswriters = Lists.newArrayList();
	public SortRealignedAssemblies(ProcessingContext processContext, AssemblyEvidenceSource source) {
		super(processContext);
		this.source = source;
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
			close();
			sort(threadpool);
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
	private void sort(ExecutorService threadpool) {
		FileSystemContext fsc = processContext.getFileSystemContext();
		List<SortCallable> tasks = Lists.newArrayList();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				tasks.add(new SortCallable(processContext,
						fsc.getAssemblyRemoteUnsortedVcfForChr(source.getFileIntermediateDirectoryBasedOn(), seq.getSequenceName()),
						fsc.getAssemblyRemoteVcfForChr(source.getFileIntermediateDirectoryBasedOn(), seq.getSequenceName()),
						VariantContextDirectedBreakpoint.ByRemoteBreakendLocationStartRaw(processContext)));
			}
		} else {
			tasks.add(new SortCallable(processContext,
					fsc.getAssemblyRemoteUnsortedVcf(source.getFileIntermediateDirectoryBasedOn()),
					fsc.getAssemblyRemoteVcf(source.getFileIntermediateDirectoryBasedOn()),
					VariantContextDirectedBreakpoint.ByRemoteBreakendLocationStartRaw(processContext)));
		}
		if (threadpool != null) {
			try {
				for (Future<Void> future : threadpool.invokeAll(tasks)) {
					// throw any exception
					future.get();
				}
			} catch (InterruptedException e) {
				log.error(e);
				throw new RuntimeException("Interrupted waiting for VCF sort", e);
			} catch (ExecutionException e) {
				throw new RuntimeException("Failed to sort VCF", e.getCause());
			}
		} else {
			for (SortCallable c : tasks) {
				c.call();
			}
		}
	}
	@Override
	public void close() {
		super.close();
		asswriters.clear();
	}
	private void writeUnsortedOutput() {
		Iterator<VariantContextDirectedBreakpoint> it = Iterators.filter(source.iterator(), VariantContextDirectedBreakpoint.class);
		while (it.hasNext()) {
			VariantContextDirectedBreakpoint bp = it.next();
			asswriters.get(bp.getBreakendSummary().referenceIndex2 % asswriters.size()).add(bp);
		}
	}
	private void createUnsortedOutputWriters() {
		FileSystemContext fsc = processContext.getFileSystemContext();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				asswriters.add(processContext.getVariantContextWriter(fsc.getAssemblyRemoteUnsortedVcfForChr(source.getFileIntermediateDirectoryBasedOn(), seq.getSequenceName())));
				toClose.add(asswriters.get(asswriters.size() - 1));
			}
		} else {
			asswriters.add(processContext.getVariantContextWriter(fsc.getAssemblyRemoteUnsortedVcf(source.getFileIntermediateDirectoryBasedOn())));
			toClose.add(asswriters.get(asswriters.size() - 1));
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
				outputs.add(fsc.getAssemblyRemoteVcfForChr(source.getFileIntermediateDirectoryBasedOn(), seq.getSequenceName()));
			}
		} else {
			outputs.add(fsc.getAssemblyRemoteVcf(source.getFileIntermediateDirectoryBasedOn()));
		}
		return outputs;
	}
	@Override
	public List<File> getTemporary() {
		List<File> files = Lists.newArrayList();
		FileSystemContext fsc = processContext.getFileSystemContext();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				files.add(fsc.getAssemblyRemoteUnsortedVcfForChr(source.getFileIntermediateDirectoryBasedOn(), seq.getSequenceName()));
			}
		} else {
			files.add(fsc.getAssemblyRemoteUnsortedVcf(source.getFileIntermediateDirectoryBasedOn()));
		}
		return files;
	}
}
