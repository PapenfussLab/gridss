package au.edu.wehi.idsv;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ExecutorService;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.util.FileHelper;
import au.edu.wehi.idsv.validation.BreakpointFilterTracker;
import au.edu.wehi.idsv.validation.OrderAssertingIterator;
import au.edu.wehi.idsv.validation.PairedEvidenceTracker;
import au.edu.wehi.idsv.vcf.VcfFileUtil;
import au.edu.wehi.idsv.visualisation.TrackedBuffer;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;


/**
 * Calls structural variant
 * @author Daniel Cameron
 *
 */
public class VariantCaller extends EvidenceProcessorBase {
	private static final Log log = Log.getInstance(VariantCaller.class);
	//private final EvidenceToCsv evidenceDump;
	public VariantCaller(ProcessingContext context, File output, List<SAMEvidenceSource> samEvidence, AssemblyEvidenceSource assemblyEvidence) {
		super(context, output, samEvidence, assemblyEvidence);
		//this.evidenceDump = evidenceDump;
	}
	@Override
	public void process(ExecutorService threadpool) {
		callBreakends(threadpool);
		annotateBreakpoints(threadpool);
	}
	public void callBreakends(ExecutorService threadpool) {
		log.info("Identifying Breakpoints");
		try {
			File outfile = processContext.getFileSystemContext().getBreakpointVcf(output);
			if (outfile.exists()) {
				log.debug("Skipping breakpoint identification: output file exists");
			} else {
				CloseableIterator<DirectedEvidence> evidenceIt = null;
				File unsorted = FileSystemContext.getWorkingFileFor(outfile);
				File sorted = FileSystemContext.getWorkingFileFor(outfile, "sorted.");
				unsorted.delete();
				sorted.delete();
				try {
					boolean assemblyOnly = processContext.getVariantCallingParameters().callOnlyAssemblies;
					evidenceIt = getAllEvidence(true, true, !assemblyOnly, !assemblyOnly, !assemblyOnly);
					evidenceIt = adjustEvidenceStream(evidenceIt);
					EvidenceClusterProcessor processor = new EvidenceClusterProcessor(processContext, evidenceIt);
					writeMaximalCliquesToVcf(
							processContext,
							processor,
							unsorted);
				} finally {
					CloserUtil.close(evidenceIt);
				}
				log.info("Sorting identified breakpoints");
				VcfFileUtil.sort(processContext, unsorted, sorted);
				FileHelper.move(sorted, outfile, true);
			}
			log.info("Identifying Breakpoints complete");
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			close();
		}
	}
	private CloseableIterator<DirectedEvidence> adjustEvidenceStream(CloseableIterator<DirectedEvidence> evidenceIt) {
		if (Defaults.SANITY_CHECK_ITERATORS) {
			evidenceIt = new PairedEvidenceTracker<DirectedEvidence>("VariantCaller", evidenceIt);
			evidenceIt = new OrderAssertingIterator<DirectedEvidence>(evidenceIt, DirectedEvidenceOrder.ByNatural);
		}
		// back to treating assemblies as independent evidence which do not affect SC or RP support counts
		// due to annotation quality being greater (in some cases over 10x) than the called quality due
		// to assembly evidence enlistment at BP evidence for a BP the evidence itself does not support
		//evidenceIt = new AutoClosingIterator<DirectedEvidence>(new OrthogonalEvidenceIterator(processContext.getLinear(), evidenceIt, getMaxWindowSize()), Lists.<Closeable>newArrayList(evidenceIt));
		return evidenceIt;
	}
	private List<Closeable> toClose = Lists.newArrayList();
	private CloseableIterator<IdsvVariantContext> getAllCalledVariants() {
		CloseableIterator<IdsvVariantContext> it = getVariants(processContext.getFileSystemContext().getBreakpointVcf(output));
		return it;
		//return new AsyncBufferedIterator<IdsvVariantContext>(it, "CalledBreakPoints");
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
	private CloseableIterator<IdsvVariantContext> getVariants(File file) {
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
		return new AutoClosingIterator<IdsvVariantContext>(idsvIt, toClose);
	}
	private static void writeMaximalCliquesToVcf(ProcessingContext processContext, Iterator<VariantContextDirectedEvidence> it, File vcf) {
		final ProgressLogger writeProgress = new ProgressLogger(log);
		VariantContextWriter vcfWriter = null;
		try {
			vcfWriter = processContext.getVariantContextWriter(vcf, false);
			log.info("Start calling maximal cliques for ", vcf);
			while (it.hasNext()) {
				VariantContextDirectedEvidence loc = it.next();
				if (loc.getPhredScaledQual() >= processContext.getVariantCallingParameters().minScore || processContext.getVariantCallingParameters().writeFiltered) {
					// If we're under min score with all possible evidence allocated, we're definitely going to fail
					// when we restrict evidence to single breakpoint support
					vcfWriter.add(loc);
				}
				writeProgress.record(processContext.getDictionary().getSequence(loc.getBreakendSummary().referenceIndex).getSequenceName(), loc.getBreakendSummary().start);
			}
			log.info("Complete calling maximal cliques for ", vcf);
		} finally {
			if (vcfWriter != null) vcfWriter.close();
		}
	}
	private int getMaxWindowSize() {
		// TODO: technically we also need to know max assembly length as the entire assembly could be a microhomology.
		int maxSize = 0;
		for (EvidenceSource source : samEvidence) {
			SAMEvidenceSource samSource = (SAMEvidenceSource)source;
			maxSize = Math.max(samSource.getMaxConcordantFragmentSize(), Math.max(samSource.getMaxReadLength(), samSource.getMaxReadMappedLength()));
		}
		maxSize = Math.max(maxSize, assemblyEvidence.getAssemblyWindowSize()) + processContext.getVariantCallingParameters().maxBreakendHomologyLength;
		return maxSize + 2 * (processContext.getVariantCallingParameters().breakendMargin + 1);
	}
	public void annotateBreakpoints(ExecutorService threadpool) {
		log.info("Annotating Calls");
		List<List<SAMEvidenceSource>> byCategory = Lists.newArrayList();
		while (byCategory.size() < processContext.getCategoryCount()) {
			byCategory.add(new ArrayList<SAMEvidenceSource>());
		}
		final int maxWindowSize = getMaxWindowSize();
		for (SAMEvidenceSource source : samEvidence) {
			byCategory.get(source.getSourceCategory()).add(source);
		}
		VariantContextWriter vcfWriter = null;
		CloseableIterator<IdsvVariantContext> it = null;
		CloseableIterator<DirectedEvidence> evidenceIt = null;
		List<ReferenceCoverageLookup> coverage = Lists.newArrayList();
		try {
			for (int i = 0; i < byCategory.size(); i++) {
				coverage.add(getReferenceLookup(byCategory.get(i), maxWindowSize));
			}
			File working = FileSystemContext.getWorkingFileFor(output);
			vcfWriter = processContext.getVariantContextWriter(working, true);
			it = getAllCalledVariants();
			Iterator<VariantContextDirectedEvidence> breakendIt = Iterators.filter(it, VariantContextDirectedEvidence.class);
			if (Defaults.SANITY_CHECK_ITERATORS) {
				breakendIt = new OrderAssertingIterator<VariantContextDirectedEvidence>(breakendIt, IdsvVariantContext.ByLocationStart);
			}
			// reorder from VCF order to breakend position order
			breakendIt = new DirectEvidenceWindowedSortingIterator<VariantContextDirectedEvidence>(processContext, maxWindowSize, breakendIt);
			processContext.registerBuffer(output.getName(), ((TrackedBuffer)breakendIt));
			if (Defaults.SANITY_CHECK_ITERATORS) {
				breakendIt = new OrderAssertingIterator<VariantContextDirectedEvidence>(breakendIt, DirectedEvidenceOrder.ByNatural);
			}
			evidenceIt = getAllEvidence(true, true, true, true, true);
			evidenceIt = adjustEvidenceStream(evidenceIt);
			breakendIt = new SequentialEvidenceAnnotator(processContext, breakendIt, evidenceIt, maxWindowSize, true, Runtime.getRuntime().availableProcessors(), threadpool);
			processContext.registerBuffer(output.getName(), ((TrackedBuffer)breakendIt));
			// breakpoint position is recalculated, so we need to resort again
			breakendIt = new DirectEvidenceWindowedSortingIterator<VariantContextDirectedEvidence>(processContext, maxWindowSize, breakendIt);
			processContext.registerBuffer(output.getName() + ".reorder", ((TrackedBuffer)breakendIt));
			if (Defaults.SANITY_CHECK_ITERATORS) {
				breakendIt = new OrderAssertingIterator<VariantContextDirectedEvidence>(breakendIt, DirectedEvidenceOrder.ByNatural);
			}
			//breakendIt = new AsyncBufferedIterator<VariantContextDirectedEvidence>(breakendIt, "Annotator-SV");
			breakendIt = new SequentialCoverageAnnotator(processContext, breakendIt, coverage);
			//breakendIt = new AsyncBufferedIterator<VariantContextDirectedEvidence>(breakendIt, "Annotator-Coverage");
			// Resort back into VCF sort order
			breakendIt = new VariantContextWindowedSortingIterator<VariantContextDirectedEvidence>(processContext, maxWindowSize, breakendIt);
			processContext.registerBuffer(output.getName() + ".vcf", ((TrackedBuffer)breakendIt));
			if (Defaults.SANITY_CHECK_ITERATORS) {
				breakendIt = new OrderAssertingIterator<VariantContextDirectedEvidence>(breakendIt, IdsvVariantContext.ByLocationStart);
				breakendIt = new BreakpointFilterTracker<VariantContextDirectedEvidence>(breakendIt, false);
			}
			while (breakendIt.hasNext()) {
				VariantContextDirectedEvidence variant = breakendIt.next();
				assert(variant.isValid());
				if (variant.isNotFiltered() || processContext.getVariantCallingParameters().writeFiltered) {
					// add confidence filter
					variant = processContext.getVariantCallingParameters().applyConfidenceFilter(processContext, variant);
					vcfWriter.add(variant);
				}
			}
			CloserUtil.close(vcfWriter);
			CloserUtil.close(breakendIt);
			vcfWriter = null;
			FileHelper.move(working, output, true);
			log.info("Variant calls written to ", output);
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			CloserUtil.close(vcfWriter);
			for (int i = 0; i < coverage.size(); i++) {
				CloserUtil.close(coverage.get(i));
			}
			CloserUtil.close(it);
			CloserUtil.close(evidenceIt);
		}
	}
	public ReferenceCoverageLookup getReferenceLookup(List<SAMEvidenceSource> input, int windowSize) {
		List<ReferenceCoverageLookup> lookup = Lists.newArrayList();
		for (SAMEvidenceSource s : input) {
			// one read-ahead thread per input file
			CloseableIterator<SAMRecord> it = new AsyncBufferedIterator<SAMRecord>(processContext.getSamReaderIterator(s.getSourceFile(), SortOrder.coordinate), s.getSourceFile().getName() + "-Coverage", gridss.Defaults.ASYNC_BUFFERS, gridss.Defaults.ASYNC_BUFFER_SIZE);
			it = new ProgressLoggingSAMRecordIterator(it, new ProgressLogger(log));
			toClose.add(it);
			SequentialReferenceCoverageLookup sourceLookup = new SequentialReferenceCoverageLookup(it, s.getMetrics().getIdsvMetrics(), s.getReadPairConcordanceCalculator(), windowSize);
			processContext.registerBuffer(s.getSourceFile().getName(), sourceLookup);
			lookup.add(sourceLookup);
		}
		return new AggregateReferenceCoverageLookup(lookup);
	}
}
