package au.edu.wehi.idsv;

import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.AutoClosingMergedIterator;
import au.edu.wehi.idsv.util.FileHelper;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;


/**
 * Calls structural variant
 * @author cameron.d
 *
 */
public class VariantCaller extends EvidenceProcessorBase {
	private static final Log log = Log.getInstance(VariantCaller.class);
	public VariantCaller(ProcessingContext context, File output, List<SAMEvidenceSource> samEvidence, AssemblyReadPairEvidenceSource assemblyEvidence) {
		super(context, output, samEvidence, assemblyEvidence);
	}
	@Override
	public void process() {
		callBreakends(null);
		annotateBreakpoints();
		writeOutput();
	}
	private static class WriteMaximalCliquesForChromosome extends EvidenceProcessorBase implements Callable<Void> {
		private int chri;
		private int chrj;
		public WriteMaximalCliquesForChromosome(ProcessingContext context, File output, List<SAMEvidenceSource> samEvidence, AssemblyReadPairEvidenceSource assemblyEvidence, int chri, int chrj) {
			super(context, output, samEvidence, assemblyEvidence);
			this.chri = chri;
			this.chrj = chrj;
		}
		@Override
		public Void call() {
			final SAMSequenceDictionary dict = processContext.getReference().getSequenceDictionary();
			final String iname = dict.getSequence(chri).getSequenceName();
			final String jname = dict.getSequence(chrj).getSequenceName();
			final String task = "Identifying Breakpoints between "  + iname + " and " + jname;
			final File outfile = processContext.getFileSystemContext().getBreakpointVcf(output, iname, jname);
			if (outfile.exists()) {
				log.debug("Skipping " + task + ": output file exists");
				return null;
			}
			CloseableIterator<DirectedEvidence> evidenceIt = null;
			try {
				FileSystemContext.getWorkingFileFor(outfile).delete();
				log.info("Start " + task);
				EvidenceClusterSubsetProcessor processor = new EvidenceClusterSubsetProcessor(processContext, chri, chrj);
				evidenceIt = getEvidenceForChr(
						processContext.getVariantCallingParameters().callOnlyAssemblies,
						false,
						true,
						true,
						true,
						chri,
						chrj); 
				writeMaximalCliquesToVcf(processContext,
						processor,
						FileSystemContext.getWorkingFileFor(outfile),
						evidenceIt);
				FileHelper.move(FileSystemContext.getWorkingFileFor(outfile), outfile, true);
				log.info("Complete " + task);
			} catch (Exception e) {
				log.error(e, "Error " + task);
				throw new RuntimeException("Error " + task, e);
			} finally {
				CloserUtil.close(evidenceIt);
			}
			return null;
		}
		@Override
		public void process() {
			call();
		}
	}
	public void callBreakends(ExecutorService threadpool) {
		log.info("Identifying Breakpoints");
		try {
			if (processContext.shouldProcessPerChromosome()) {
				List<WriteMaximalCliquesForChromosome> workers = new ArrayList<>();
				final SAMSequenceDictionary dict = processContext.getReference().getSequenceDictionary();
				for (int i = 0; i < dict.size(); i++) {					
					for (int j = i; j < dict.size(); j++) {
						workers.add(new WriteMaximalCliquesForChromosome(processContext, output, samEvidence, assemblyEvidence, i, j));
					}
				}
				if (threadpool != null) {
					log.info("Identifying Breakpoints in parallel");
					try {
						for (Future<Void> future : threadpool.invokeAll(workers)) {
							// throw exception from worker thread here
							future.get();
						}
					} catch (InterruptedException e) {
						log.error(e, "Interrupted Identifying Breakpoints");
						throw new RuntimeException(e);
					} catch (ExecutionException e) {
						log.error(e, "Error Identifying Breakpoints");
						throw new RuntimeException(e);
					}
				} else {
					log.info("Identifying Breakpoints sequentially");
					for (WriteMaximalCliquesForChromosome c : workers) {
						c.call();
					}
				}
			} else {
				File outfile = processContext.getFileSystemContext().getBreakpointVcf(output);
				if (outfile.exists()) {
					log.debug("Skipping breakpoint identification: output file exists");
				} else {
					CloseableIterator<DirectedEvidence> evidenceIt = null;
					try {
						EvidenceClusterProcessor processor = new EvidenceClusterProcessor(processContext);
						// TODO: linear pass over breakpoints by including all remote evidence
						evidenceIt = getAllEvidence(processContext.getVariantCallingParameters().callOnlyAssemblies,
								false,
								true,
								true,
								true);
						writeMaximalCliquesToVcf(
								processContext,
								processor,
								FileSystemContext.getWorkingFileFor(outfile),
								evidenceIt);
						FileHelper.move(FileSystemContext.getWorkingFileFor(outfile), outfile, true);
					} finally {
						CloserUtil.close(evidenceIt);
					}
				}
			}
			log.info("Identifying Breakpoints complete");
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			close();
		}
	}
	private List<Closeable> toClose = new ArrayList<>();
	private CloseableIterator<IdsvVariantContext> getAllCalledVariants() {
		List<Iterator<IdsvVariantContext>> variants = new ArrayList<>();
		if (processContext.shouldProcessPerChromosome()) {
			SAMSequenceDictionary dict = processContext.getReference().getSequenceDictionary();
			for (int i = 0; i < dict.size(); i++) {
				final String chri = dict.getSequence(i).getSequenceName();
				for (int j = i; j < dict.size(); j++) {
					final String chrj = dict.getSequence(j).getSequenceName();
					variants.add(getVariants(processContext.getFileSystemContext().getBreakpointVcf(output, chri, chrj)));
				}
			}
		} else {
			variants.add(getVariants(processContext.getFileSystemContext().getBreakpointVcf(output)));
		}
		CloseableIterator<IdsvVariantContext> merged = new AutoClosingMergedIterator<IdsvVariantContext>(variants, IdsvVariantContext.ByLocationStart);
		return new AsyncBufferedIterator<>(merged, "all breakpoints");
	}
	public void close() {
		super.close();
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
	private Iterator<IdsvVariantContext> getVariants(File file) {
		VCFFileReader vcfReader = new VCFFileReader(file, false);
		toClose.add(vcfReader);
		CloseableIterator<VariantContext> it = vcfReader.iterator();
		toClose.add(it);
		Iterator<IdsvVariantContext> idsvIt = Iterators.transform(it, new Function<VariantContext, IdsvVariantContext>() {
			@Override
			public IdsvVariantContext apply(VariantContext arg) {
				return IdsvVariantContext.create(processContext, null, arg);
			}
		});
		return idsvIt;
	}
	private static void writeMaximalCliquesToVcf(ProcessingContext processContext, EvidenceClusterProcessor processor, File vcf, Iterator<DirectedEvidence> evidenceIt) {
		final ProgressLogger writeProgress = new ProgressLogger(log);
		log.info("Loading minimal evidence set for ", vcf, " into memory.");
		while (evidenceIt.hasNext()) {
			processor.addEvidence(evidenceIt.next());
		}
		VariantContextWriter vcfWriter = null;
		try {
			vcfWriter = processContext.getVariantContextWriter(vcf, true);
			log.info("Start calling maximal cliques for ", vcf);
			Iterator<VariantContextDirectedEvidence> it = processor.iterator();
			while (it.hasNext()) {
				VariantContextDirectedEvidence loc = it.next();
				vcfWriter.add(loc);
				writeProgress.record(processContext.getDictionary().getSequence(loc.getBreakendSummary().referenceIndex).getSequenceName(), loc.getBreakendSummary().start);
			}
			log.info("Complete calling maximal cliques for ", vcf);
		} finally {
			if (vcfWriter != null) vcfWriter.close();
		}
	}
	public void annotateBreakpoints() {
		annotateBreakpoints(null);
	}
	public void annotateBreakpoints(BreakendAnnotator annotator) {
		log.info("Annotating Calls");
		List<SAMEvidenceSource> normal = new ArrayList<>();
		List<SAMEvidenceSource> tumour = new ArrayList<>();
		int maxFragmentSize = assemblyEvidence.getMaxConcordantFragmentSize();
		for (EvidenceSource source : samEvidence) {
			maxFragmentSize = Math.max(source.getMaxConcordantFragmentSize(), maxFragmentSize);
			SAMEvidenceSource samSource = (SAMEvidenceSource)source;
			if (samSource.isTumour()) {
				tumour.add(samSource);
			} else {
				normal.add(samSource);
			}
		}
		VariantContextWriter vcfWriter = null;
		ReferenceCoverageLookup normalCoverage = null;
		ReferenceCoverageLookup tumourCoverage = null;
		CloseableIterator<IdsvVariantContext> it = null;
		CloseableIterator<DirectedEvidence> evidenceIt = null;
		try {
			vcfWriter = processContext.getVariantContextWriter(FileSystemContext.getWorkingFileFor(output), true);
			normalCoverage = getReferenceLookup(normal, 4 * maxFragmentSize);
			tumourCoverage = getReferenceLookup(tumour, 4 * maxFragmentSize);
			BreakendAnnotator referenceAnnotator  = new SequentialCoverageAnnotator(processContext, normalCoverage, tumourCoverage);
			evidenceIt = getAllEvidence(true, true, true, true, true);
			BreakendAnnotator evidenceAnnotator = new SequentialEvidenceAnnotator(processContext, evidenceIt);
			it = getAllCalledVariants();
			while (it.hasNext()) {
				IdsvVariantContext rawVariant = it.next();
				if (rawVariant instanceof VariantContextDirectedEvidence) {
					VariantContextDirectedEvidence annotatedVariant = (VariantContextDirectedEvidence)rawVariant;
					if (rawVariant.isValid()) {
						annotatedVariant = referenceAnnotator.annotate((VariantContextDirectedEvidence)rawVariant);
						annotatedVariant = evidenceAnnotator.annotate(annotatedVariant);
						if (annotator != null) {
							annotatedVariant = annotator.annotate(annotatedVariant);
						}
						if (!annotatedVariant.isFiltered() || Defaults.WRITE_FILTERED_CALLS) {
							vcfWriter.add(annotatedVariant);
						}
					}
				}
			}
			CloserUtil.close(vcfWriter);
			FileHelper.move(FileSystemContext.getWorkingFileFor(output), output, true);
			log.info("Variant calls written to ", output);
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			CloserUtil.close(vcfWriter);
			CloserUtil.close(normalCoverage);
			CloserUtil.close(tumourCoverage);
			CloserUtil.close(it);
			CloserUtil.close(evidenceIt);
		}
	}
	public ReferenceCoverageLookup getReferenceLookup(List<SAMEvidenceSource> input, int windowSize) {
		List<ReferenceCoverageLookup> lookup = new ArrayList<>();
		for (SAMEvidenceSource s : input) {
			// one read-ahead thread per input file
			CloseableIterator<SAMRecord> it = new AsyncBufferedIterator<>(processContext.getSamReaderIterator(s.getSourceFile(), SortOrder.coordinate), s.getSourceFile().getName() + " reference coverage");
			toClose.add(it);
			lookup.add(new SequentialReferenceCoverageLookup(it, s.getReadPairConcordanceCalculator(), windowSize));
		}
		return new AggregateReferenceCoverageLookup(lookup);
	}
	private void writeOutput() {
		log.info("Outputting calls to ", output);
	}
}
