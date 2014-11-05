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

import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.concurrent.Callable;

import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.IntermediateFileUtil;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.util.FileHelper;

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
	/**
	 * VCF Sort task
	 * @author cameron.d
	 *
	 */
	public static class SortCallable implements Callable<Void> {
		private final ProcessingContext processContext;
		private final File input;
		private final File output;
		private final Comparator<VariantContext> sortComparator;
		public SortCallable(ProcessingContext processContext, File input, File output, Comparator<VariantContext> sortComparator) {
			this.processContext = processContext;
			this.input = input;
			this.output = output;
			this.sortComparator = sortComparator;
		}
		@Override
		public Void call() throws IOException {
			if (IntermediateFileUtil.checkIntermediate(output, input)) {
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
				reader.close();
				collection.doneAdding();
				writer = processContext.getVariantContextWriter(FileSystemContext.getWorkingFileFor(output));
		    	wit = collection.iterator();
				while (wit.hasNext()) {
					writer.add(wit.next());
				}
				wit.close();
				writer.close();
				collection.cleanup();
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
