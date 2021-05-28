package au.edu.wehi.idsv;

import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.FileHelper;
import au.edu.wehi.idsv.vcf.VcfFileUtil;
import com.google.common.util.concurrent.MoreExecutors;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;


/**
 * Calls structural variant
 * @author Daniel Cameron
 *
 */
public class VariantCaller {
	private static final Log log = Log.getInstance(VariantCaller.class);
	private final ProcessingContext processContext;
	private final List<SAMEvidenceSource> samEvidence;
	private final List<AssemblyEvidenceSource> assemblyEvidence;
	//private final EvidenceToCsv evidenceDump;
	public VariantCaller(ProcessingContext context, List<SAMEvidenceSource> samEvidence, List<AssemblyEvidenceSource> assemblyEvidence) {
		this.processContext = context;
		this.samEvidence = samEvidence;
		this.assemblyEvidence = assemblyEvidence;
	}
	public void callBreakends(File vcf, ExecutorService threadpool) throws IOException {
		samEvidence.stream().forEach(ses -> ses.assertPreprocessingComplete());
		for (AssemblyEvidenceSource aes : assemblyEvidence) {
			aes.assertPreprocessingComplete();
		}
		log.info("Identifying Breakpoints");
		if (threadpool == null) {
			threadpool = MoreExecutors.newDirectExecutorService();
		}
		AggregateEvidenceSource es = new AggregateEvidenceSource(processContext, samEvidence, assemblyEvidence, SAMEvidenceSource.EvidenceSortOrder.EvidenceStartPosition);
		List<QueryInterval[]> chunks = processContext.getReference().getIntervals(processContext.getConfig().chunkSize, processContext.getConfig().chunkSequenceChangePenalty);
		List<File> calledChunk = new ArrayList<>();
		List<Future<Void>> tasks = new ArrayList<>();
		
		for (int i = 0; i < chunks.size(); i++) {
			QueryInterval[] chunk = chunks.get(i);
			File f = processContext.getFileSystemContext().getVariantCallChunkVcf(vcf, i);
			int chunkNumber = i;
			calledChunk.add(f);
			if (!f.exists()) {
				tasks.add(threadpool.submit(() -> { callChunk(f, es, chunkNumber, chunk); return null; }));
			}
		}
		runTasks(tasks);
		
		log.info("Merging identified breakpoints");
		File mergedOut = FileSystemContext.getWorkingFileFor(vcf, "gridss.merged.");
		VcfFileUtil.concat(processContext.getReference().getSequenceDictionary(), calledChunk, mergedOut);
		
		log.info("Sorting identified breakpoints");
		VcfFileUtil.sort(processContext, mergedOut, vcf);
		// clean up chunked
		if (gridss.Defaults.DELETE_TEMPORARY_FILES) {
			for (File f : calledChunk) {
				FileHelper.delete(f, true);
			}
			FileHelper.delete(mergedOut, true);
		}
	}
	private void runTasks(List<Future<Void>> tasks) {
		// Run as many tasks as we can before dying
		Exception firstException = null;
		for (Future<Void> f : tasks) {
			try {
				f.get();
			} catch (Exception e) {
				if (firstException == null) {
					firstException = e;
				}
			}
		}
		if (firstException != null) {
			log.error(firstException, "Fatal error during breakpoint identification ");
			throw new RuntimeException(firstException);
		}
	}
	private void callChunk(File output, AggregateEvidenceSource es, int chunkNumber, QueryInterval[] chunk) {
		try {
			String chunkMsg = String.format("chunk %d (%s:%d-%s:%d)", chunkNumber,
					processContext.getDictionary().getSequence(chunk[0].referenceIndex).getSequenceName(), chunk[0].start,
					processContext.getDictionary().getSequence(chunk[chunk.length - 1].referenceIndex).getSequenceName(), chunk[chunk.length - 1].end);
			String msg = "calling maximal cliques in " + chunkMsg;
			File tmp = new File(output.getParent(), "gridss.tmp." + output.getName());
			try (VariantCallIterator rawit = new VariantCallIterator(es, chunk, chunkNumber)) {
				try (VariantContextWriter vcfWriter = processContext.getVariantContextWriter(tmp, false)) {
					log.info("Start ", msg);
					try (AsyncBufferedIterator<VariantContextDirectedEvidence> it = new AsyncBufferedIterator<>(rawit, "VariantCaller " + chunkMsg)) {
						while (it.hasNext()) {
							VariantContextDirectedEvidence loc = it.next();
							boolean hardFiltered = processContext.getVariantCallingParameters().isHardFilteredBeforeAnnotation(loc);
							if (!hardFiltered || processContext.getVariantCallingParameters().writeFiltered) {
								vcfWriter.add(loc);
							}
						}
					}
				}
			}
			try {
				FileHelper.move(tmp, output, true);
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
			log.info("Complete ", msg);
			if (gridss.Defaults.DEFENSIVE_GC) {
				log.info("Requesting defensive GC to ensure OS file handles are closed");
				System.gc();
				System.runFinalization();
			}
		} catch (OutOfMemoryError oom) {
			log.error(oom);
			System.exit(1);
		}
	}
}
