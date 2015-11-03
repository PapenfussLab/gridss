package au.edu.wehi.idsv.vcf;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFRecordCodec;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.concurrent.Callable;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.IdsvVariantContext;
import au.edu.wehi.idsv.IntermediateFileUtil;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.util.FileHelper;
import au.edu.wehi.idsv.validation.OrderAssertingIterator;

import com.google.common.collect.ImmutableList;

public class VcfFileUtil {
	private static final Log log = Log.getInstance(VcfFileUtil.class);
	/**
	 * Sorts a VCF according to the given sort order
	 * @param input unsorted input
	 * @param output sorted output to write
	 * @param sortComparator sort order
	 * @throws IOException 
	 */
	public static void sort(ProcessingContext processContext, File input, File output, Comparator<VariantContext> sortComparator) {
		try {
			new SortCallable(processContext, input, output, sortComparator).call();
		} catch (IOException e) {
			log.error(log);
			throw new RuntimeException(e);
		}
	}
	public static void sort(ProcessingContext processContext, File input, File output) {
		try {
			new SortCallable(processContext, input, output).call();
		} catch (IOException e) {
			log.error(log);
			throw new RuntimeException(e);
		}
	}
	/**
	 * VCF Sort task
	 * @author Daniel Cameron
	 *
	 */
	public static class SortCallable implements Callable<Void> {
		private final ProcessingContext processContext;
		private final File input;
		private final File output;
		private final Comparator<VariantContext> sortComparator;
		private final boolean indexed;
		public SortCallable(ProcessingContext processContext, File input, File output) {
			this(processContext, input, output, IdsvVariantContext.VariantContextByLocationStart(processContext.getDictionary()));
		}
		public SortCallable(ProcessingContext processContext, File input, File output, Comparator<VariantContext> sortComparator) {
			this(processContext, input, output, sortComparator, false);
		}
		private SortCallable(ProcessingContext processContext, File input, File output, Comparator<VariantContext> sortComparator, boolean writeIndex) {
			this.processContext = processContext;
			this.input = input;
			this.output = output;
			this.sortComparator = sortComparator;
			this.indexed = writeIndex;
		}
		@Override
		public Void call() throws IOException {
			if (IntermediateFileUtil.checkIntermediate(output, input, processContext.getConfig().ignoreFileTimestamps)) {
				log.info("Not sorting as output already exists: " + output);
				return null;
			}
			log.info("Sorting " + input);
			VCFFileReader reader = null;
			VariantContextWriter writer = null;
			CloseableIterator<VariantContext> rit = null;
			CloseableIterator<VariantContext> wit = null;
			SortingCollection<VariantContext> collection = null;
			if (FileSystemContext.getWorkingFileFor(output).exists()) FileSystemContext.getWorkingFileFor(output).delete();
			try {
				reader = new VCFFileReader(input, false);
				VCFHeader header = reader.getFileHeader();
				rit = reader.iterator();
				collection = SortingCollection.newInstance(
						VariantContext.class,
						new VCFRecordCodec(header),
						sortComparator,
						processContext.getFileSystemContext().getMaxBufferedRecordsPerFile(),
						processContext.getFileSystemContext().getTemporaryDirectory());
				while (rit.hasNext()) {
					collection.add(rit.next());
				}
				rit.close();
				rit = null;
				reader.close();
				reader = null;
				collection.doneAdding();
				writer = processContext.getVariantContextWriter(FileSystemContext.getWorkingFileFor(output), indexed);
		    	wit = collection.iterator();
		    	if (Defaults.SANITY_CHECK_ITERATORS) {
					wit = new AutoClosingIterator<VariantContext>(new OrderAssertingIterator<VariantContext>(wit, sortComparator), ImmutableList.<Closeable>of(wit));
				}
				while (wit.hasNext()) {
					writer.add(wit.next());
				}
				wit.close();
				wit = null;
				writer.close();
				writer = null;
				collection.cleanup();
				collection = null;
				FileHelper.move(FileSystemContext.getWorkingFileFor(output), output, true);
			} finally {
				CloserUtil.close(writer);
				CloserUtil.close(wit);
				CloserUtil.close(rit);
				CloserUtil.close(reader);
				if (collection != null) collection.cleanup();
				if (FileSystemContext.getWorkingFileFor(output).exists()) FileSystemContext.getWorkingFileFor(output).delete();
			}
			return null;
		}
	}
}
