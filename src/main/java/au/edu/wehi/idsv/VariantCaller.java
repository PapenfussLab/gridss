package au.edu.wehi.idsv;

import java.io.File;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ExecutorService;

import au.edu.wehi.idsv.validation.OrderAssertingIterator;
import au.edu.wehi.idsv.validation.PairedEvidenceTracker;
import au.edu.wehi.idsv.vcf.VcfFileUtil;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;


/**
 * Calls structural variant
 * @author Daniel Cameron
 *
 */
public class VariantCaller {
	private static final Log log = Log.getInstance(VariantCaller.class);
	private final ProcessingContext processContext;
	private final List<SAMEvidenceSource> samEvidence;
	private final AssemblyEvidenceSource assemblyEvidence;
	//private final EvidenceToCsv evidenceDump;
	public VariantCaller(ProcessingContext context, List<SAMEvidenceSource> samEvidence, AssemblyEvidenceSource assemblyEvidence) {
		this.processContext = context;
		this.samEvidence = samEvidence;
		this.assemblyEvidence = assemblyEvidence;
	}
	public void callBreakends(File vcf, ExecutorService threadpool) {
		log.info("Identifying Breakpoints");
		File tmp = FileSystemContext.getWorkingFileFor(vcf);
		CloseableIterator<DirectedEvidence> evidenceIt = null;
		try {
			AggregateEvidenceSource es = new AggregateEvidenceSource(
					processContext,
					processContext.getVariantCallingParameters().callOnlyAssemblies ? Collections.emptyList() : samEvidence, assemblyEvidence);
			evidenceIt = adjustEvidenceStream(evidenceIt);
			EvidenceClusterProcessor processor = new EvidenceClusterProcessor(threadpool, es);
			writeMaximalCliquesToVcf(
					processContext,
					processor,
					tmp);
		} finally {
			CloserUtil.close(evidenceIt);
		}
		log.info("Sorting identified breakpoints");
		VcfFileUtil.sort(processContext, tmp, vcf);
	}
	private CloseableIterator<DirectedEvidence> adjustEvidenceStream(CloseableIterator<DirectedEvidence> evidenceIt) {
		if (Defaults.SANITY_CHECK_ITERATORS) {
			evidenceIt = new PairedEvidenceTracker<DirectedEvidence>("VariantCaller", evidenceIt);
			evidenceIt = new OrderAssertingIterator<DirectedEvidence>(evidenceIt, DirectedEvidenceOrder.ByNatural);
		}
		return evidenceIt;
	}
	private static void writeMaximalCliquesToVcf(ProcessingContext processContext, Iterator<VariantContextDirectedEvidence> it, File vcf) {
		final ProgressLogger writeProgress = new ProgressLogger(log);
		try (VariantContextWriter vcfWriter = processContext.getVariantContextWriter(vcf, false)) {
			log.info("Start calling maximal cliques for ", vcf);
			while (it.hasNext()) {
				VariantContextDirectedEvidence loc = it.next();
				if (loc.getBreakendQual() >= processContext.getVariantCallingParameters().minScore || processContext.getVariantCallingParameters().writeFiltered) {
					// If we're under min score with all possible evidence allocated, we're definitely going to fail
					// when we restrict evidence to single breakpoint support
					vcfWriter.add(loc);
				}
				writeProgress.record(processContext.getDictionary().getSequence(loc.getBreakendSummary().referenceIndex).getSequenceName(), loc.getBreakendSummary().start);
			}
			log.info("Complete calling maximal cliques for ", vcf);
		}
	}
}
