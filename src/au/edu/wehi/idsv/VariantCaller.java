package au.edu.wehi.idsv;

import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;


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
		Iterator<IdsvVariantContext> merged = Iterators.mergeSorted(variants, IdsvVariantContext.ByLocationStart);
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
		log.info("Loading minimal evidence set for ", vcf, " into memory.");
		while (evidenceIt.hasNext()) {
			processor.addEvidence(evidenceIt.next());
		}
		VariantContextWriter vcfWriter = null;
		try {
			vcfWriter = processContext.getVariantContextWriter(vcf);
			toClose.add(vcfWriter);
			log.info("Calling maximal cliques for ", vcf);
			Iterator<VariantContextDirectedEvidence> it = processor.iterator();
			while (it.hasNext()) {
				VariantContextDirectedEvidence loc = it.next();
				vcfWriter.add(loc);
				writeProgress.record(processContext.getDictionary().getSequence(loc.getBreakendSummary().referenceIndex).getSequenceName(), loc.getBreakendSummary().start);
			}
		} finally {
			if (vcfWriter != null) vcfWriter.close();
		}
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
		try {
			final VariantContextWriter vcfWriter = processContext.getVariantContextWriter(output);
			toClose.add(vcfWriter);
			final SequentialReferenceCoverageLookup normalCoverage = getReferenceLookup(normalFiles);
			toClose.add(normalCoverage);
			final SequentialReferenceCoverageLookup tumourCoverage = getReferenceLookup(tumourFiles);
			toClose.add(tumourCoverage);
			SequentialBreakendAnnotator annotator = new SequentialBreakendAnnotator(processContext, normalCoverage, tumourCoverage, getAllEvidence());
			Iterator<IdsvVariantContext> it = getAllCalledVariants();
			while (it.hasNext()) {
				IdsvVariantContext rawVariant = it.next();
				if (rawVariant instanceof VariantContextDirectedEvidence && ((VariantContextDirectedEvidence)rawVariant).isValid()) {
					VariantContextDirectedEvidence annotatedVariant = annotator.annotate((VariantContextDirectedEvidence)rawVariant);
					vcfWriter.add(annotatedVariant);
				}
			}
			vcfWriter.close();
			normalCoverage.close();
			tumourCoverage.close();
			log.info("Variant calls written to ", output);
		} finally {
			close();
		}
	}
	public SequentialReferenceCoverageLookup getReferenceLookup(List<File> samFiles) {
		List<CloseableIterator<SAMRecord>> toMerge = Lists.newArrayList();
		for (File f : samFiles) {
			SamReader reader = processContext.getSamReader(f);
			toClose.add(reader);
			CloseableIterator<SAMRecord> it = processContext.getSamReaderIterator(reader, SortOrder.coordinate);
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