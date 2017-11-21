package gridss.cmdline;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;

import org.broadinstitute.barclay.argparser.Argument;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;

import au.edu.wehi.idsv.DirectEvidenceWindowedSortingIterator;
import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.IdsvVariantContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.VariantContextWindowedSortingIterator;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.util.FileHelper;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import picard.cmdline.StandardOptionDefinitions;

/**
 * Base class used to transform a VCF breakpoint call set given the full evidence available.
 * 
 * 
 * @author Daniel Cameron
 *
 */
public abstract class VcfTransformCommandLineProgram extends FullEvidenceCommandLineProgram {
	private static final Log log = Log.getInstance(VcfTransformCommandLineProgram.class);
	@Argument(shortName="VCF", doc="Coordinate sorted VCF file")
    public File INPUT_VCF;
	@Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="VCF structural variation calls.")
    public File OUTPUT_VCF;
	public abstract CloseableIterator<VariantContextDirectedBreakpoint> iterator(CloseableIterator<VariantContextDirectedBreakpoint> calls, ExecutorService threadpool);
	@Override
	public int doWork(ExecutorService threadpool) throws IOException, InterruptedException, ExecutionException {
		IOUtil.assertFileIsReadable(ASSEMBLY);
		IOUtil.assertFileIsReadable(INPUT_VCF);
		IOUtil.assertFileIsWritable(OUTPUT_VCF);
		log.info("Annotating variants in " + INPUT_VCF);
		try (CloseableIterator<VariantContextDirectedBreakpoint> it = iterator(getBreakpoints(INPUT_VCF), threadpool)) {
			saveVcf(OUTPUT_VCF, getAllCalls(INPUT_VCF, it));
		}
		log.info("Annotated variants written to " + OUTPUT_VCF);
		return 0;
	}
	public CloseableIterator<VariantContextDirectedBreakpoint> getBreakpoints(File file) {
		VCFFileReader vcfReader = new VCFFileReader(file, false);
		CloseableIterator<VariantContext> it = vcfReader.iterator();
		Iterator<IdsvVariantContext> idsvIt = Iterators.transform(it, variant -> IdsvVariantContext.create(getContext(), null, variant));
		Iterator<VariantContextDirectedBreakpoint> bpit = Iterators.filter(idsvIt, VariantContextDirectedBreakpoint.class);
		// resort by evidence start
		bpit = new DirectEvidenceWindowedSortingIterator<>(getContext(), SAMEvidenceSource.maximumWindowSize(getContext(), getSamEvidenceSources(), getAssemblySource()), bpit);
		return new AutoClosingIterator<VariantContextDirectedBreakpoint>(bpit, it, vcfReader);
	}
	public Iterator<IdsvVariantContext> getAllCalls(File file, CloseableIterator<VariantContextDirectedBreakpoint> breakpointCalls) {
		VCFFileReader vcfReader = new VCFFileReader(file, false);
		CloseableIterator<VariantContext> it = vcfReader.iterator();
		Iterator<IdsvVariantContext> idsvIt = Iterators.transform(it, variant -> IdsvVariantContext.create(getContext(), null, variant));
		Iterator<IdsvVariantContext> nonbpIt = Iterators.filter(idsvIt, variant -> !(variant instanceof VariantContextDirectedBreakpoint));
		// sort back to nominal VCF position
		Iterator<VariantContextDirectedBreakpoint> bpit = new VariantContextWindowedSortingIterator<>(getContext(), SAMEvidenceSource.maximumWindowSize(getContext(), getSamEvidenceSources(), getAssemblySource()), breakpointCalls);
		Iterator<IdsvVariantContext> mergedIt = Iterators.mergeSorted(ImmutableList.of(bpit, nonbpIt), IdsvVariantContext.ByLocationStart); 
		return new AutoClosingIterator<>(mergedIt, vcfReader, it);
	}
	protected void saveVcf(File file, Iterator<IdsvVariantContext> calls) throws IOException {
		File tmp = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(file) : file;
		final ProgressLogger writeProgress = new ProgressLogger(log);
		try (VariantContextWriter vcfWriter = getContext().getVariantContextWriter(tmp, true)) {
			while (calls.hasNext()) {
				IdsvVariantContext record = calls.next();
				vcfWriter.add(record);
				writeProgress.record(record.getContig(), record.getStart());
			}
		}
		if (tmp != file) {
			FileHelper.move(tmp, file, true);
		}
	}
}
