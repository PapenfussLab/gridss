package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.FilteringIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VCFWriterUnitTest;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;

import au.edu.wehi.idsv.metrics.RelevantMetrics;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.UnmodifiableIterator;


/**
 * Calls structural variant
 * @author cameron.d
 *
 */
public class VariantCaller extends EvidenceProcessorBase {
	private static final Log log = Log.getInstance(VariantCaller.class);
	public VariantCaller(ProcessingContext context, File output, List<EvidenceSource> evidence) {
		super(context, output, evidence);
	}
	@Override
	public void process() {
		callBreakends();
		annotateBreakpoints();
		writeOutput();
	}
	public void callBreakends() {
		log.info("Identifying Breakpoints");
		try {
			if (processContext.shouldProcessPerChromosome()) {
				SAMSequenceDictionary dict = processContext.getReference().getSequenceDictionary();
				for (int i = 0; i < dict.size(); i++) {
					String chri = dict.getSequence(i).getSequenceName();
					for (int j = i; j < dict.size(); j++) {
						String chrj = dict.getSequence(j).getSequenceName();
						EvidenceClusterSubsetProcessor processor = new EvidenceClusterSubsetProcessor(processContext, i, j);
						writeMaximalCliquesToVcf(
								processor,
								processContext.getFileSystemContext().getBreakpointVcf(output, chri, chrj),
								getEvidenceForChr(i, j));
					}
				}
			} else {
				EvidenceClusterProcessor processor = new EvidenceClusterProcessor(processContext);
				writeMaximalCliquesToVcf(
						processor,
						processContext.getFileSystemContext().getBreakpointVcf(output),
						getAllEvidence());
			}
		} finally {
			close();
		}
	}
	private List<Closeable> toClose = Lists.newArrayList();
	private Iterator<IdsvVariantContext> getAllCalledVariants() {
		List<Iterator<IdsvVariantContext>> variants = Lists.newArrayList();
		try {
			if (processContext.shouldProcessPerChromosome()) {
				SAMSequenceDictionary dict = processContext.getReference().getSequenceDictionary();
				for (int i = 0; i < dict.size(); i++) {
					String chri = dict.getSequence(i).getSequenceName();
					for (int j = i; j < dict.size(); j++) {
						String chrj = dict.getSequence(j).getSequenceName();
						variants.add(getVariants(processContext.getFileSystemContext().getBreakpointVcf(output, chri, chrj)));
					}
				}
			} else {
				variants.add(getVariants(processContext.getFileSystemContext().getBreakpointVcf(output)));
			}
		} finally {
			close();
		}
		Iterator<IdsvVariantContext> merged = Iterators.mergeSorted(variants, IdsvVariantContext.ByLocation);
		return merged;
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
		VCFFileReader vcfReader = new VCFFileReader(file);
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
	private void writeMaximalCliquesToVcf(EvidenceClusterProcessor processor, File vcf, Iterator<DirectedEvidence> evidenceIt) {
		final ProgressLogger writeProgress = new ProgressLogger(log);
		log.debug("Loading minimal evidence set for ", vcf, " into memory.");
		while (evidenceIt.hasNext()) {
			processor.addEvidence(evidenceIt.next());
		}
		try (VariantContextWriter vcfWriter = processContext.getVariantContextWriter(vcf)) {
			log.debug("Calling maximal cliques for ", vcf);
			// Write out both sides of the breakend in order
			// since the first breakend is always the lower genomic coordinate
			// this will result in in-order output
			PriorityQueue<BreakpointSummary> highEnd = new PriorityQueue<BreakpointSummary>(1024, BreakendSummary.ByStartEnd);
			processor = new EvidenceClusterProcessor(processContext);
			Iterator<BreakpointSummary> it = processor.iterator();
			while (it.hasNext()) {
				BreakendSummary loc = it.next();
				if (loc instanceof BreakpointSummary) {
					// Add the remote side of the call
					highEnd.add(((BreakpointSummary)loc).remoteBreakpoint());
				}
				while (!highEnd.isEmpty() && BreakendSummary.ByStartEnd.compare(highEnd.peek(), loc) < 0) {
					// write remote calls that are before here
					writeCall(processContext, vcfWriter, highEnd.poll(), writeProgress);
				}
				// write call
				writeCall(processContext, vcfWriter, loc, writeProgress);
			}
		}
	}
	private void writeCall(final ProcessingContext processContext, final VariantContextWriter vcfWriter, BreakendSummary loc, ProgressLogger progress) {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(processContext, null)
			.breakend(loc, null) // Issue: we've lost all our extended info including the untemplated sequence
			.evidence(loc.evidence);
		VariantContext vc = builder.make(); 
		vcfWriter.add(vc);
		progress.record(vc.getChr(), vc.getStart());
	}
	public void annotateBreakpoints() {
		log.info("Annotating Calls");
		List<File> normalFiles = Lists.newArrayList();
		List<File> tumourFiles = Lists.newArrayList();
		for (EvidenceSource source : evidence) {
			if (source instanceof SAMEvidenceSource) {
				SAMEvidenceSource samSource = (SAMEvidenceSource)source;
				if (samSource.isTumour()) {
					tumourFiles.add(samSource.getSourceFile());
				} else {
					normalFiles.add(samSource.getSourceFile());
				}
			}
		}
		try (final VariantContextWriter vcfWriter = processContext.getVariantContextWriter(output)) {
			final SequentialReferenceCoverageLookup normalCoverage = getReferenceLookup(normalFiles);
			final SequentialReferenceCoverageLookup tumourCoverage = getReferenceLookup(tumourFiles);
			SequentialBreakendAnnotator annotator = new SequentialBreakendAnnotator(processContext, normalCoverage, tumourCoverage, getAllEvidence());
			Iterator<IdsvVariantContext> it = getAllCalledVariants();
			while (it.hasNext()) {
				IdsvVariantContext rawvariant = it.next();
				if (rawvariant instanceof VariantContextDirectedBreakpoint && ((VariantContextDirectedBreakpoint)rawvariant).isValid()) {
					VariantContextDirectedBreakpoint variant = (VariantContextDirectedBreakpoint)rawvariant;
					annotator.annotate(variant);
				}
			}
		} finally {
			close();
		}
	}
	public SequentialReferenceCoverageLookup getReferenceLookup(List<File> samFiles) {
		List<CloseableIterator<SAMRecord>> toMerge = Lists.newArrayList();
		for (File f : samFiles) {
			SamReader reader = processContext.getSamReaderFactory().open(f);
			toClose.add(reader);
			CloseableIterator<SAMRecord> it = reader.iterator().assertSorted(SortOrder.coordinate);
			toClose.add(it);
			it = processContext.applyCommonSAMRecordFilters(it);
			toMerge.add(it);
		}
		Iterator<SAMRecord> merged = Iterators.mergeSorted(toMerge, new SAMRecordQueryNameComparator()); 
		return new SequentialReferenceCoverageLookup(merged, 1024);
	}
	private void writeOutput() {
		log.info("Outputting calls to ", output);
	}
}
