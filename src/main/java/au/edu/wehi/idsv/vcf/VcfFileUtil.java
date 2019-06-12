package au.edu.wehi.idsv.vcf;

import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.IdsvVariantContext;
import au.edu.wehi.idsv.IntermediateFileUtil;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.FileHelper;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFRecordCodec;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.Callable;

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
			if (IntermediateFileUtil.checkIntermediate(output)) {
				log.info("Not sorting as output already exists: " + output);
				return null;
			}
			log.info("Sorting to " + output);
			SortingCollection<VariantContext> collection = null;
			File tmpout = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(output, "gridss.tmp.sorting.") : output;
			if (tmpout != output && tmpout.exists()) {
				FileHelper.delete(tmpout, true);
			}
			try {
				try (VCFFileReader reader = new VCFFileReader(input, false)) {
					VCFHeader header = reader.getFileHeader();
					try (CloseableIterator<VariantContext> rit = reader.iterator()) {
						collection = SortingCollection.newInstance(
								VariantContext.class,
								new VCFRecordCodec(header),
								sortComparator,
								processContext.getFileSystemContext().getMaxBufferedRecordsPerFile(),
								processContext.getFileSystemContext().getTemporaryDirectory().toPath());
						while (rit.hasNext()) {
							collection.add(rit.next());
						}
					}
				}
				collection.doneAdding();
				try (VariantContextWriter writer = processContext.getVariantContextWriter(tmpout, indexed)) {
			    	try (CloseableIterator<VariantContext> wit = collection.iterator()) {
						while (wit.hasNext()) {
							writer.add(wit.next());
						}
			    	}
				}
				collection.cleanup();
				collection = null;
				if (tmpout != output) {
					FileHelper.move(tmpout, output, true);
				}
			} finally {
				if (collection != null) collection.cleanup();
				if (tmpout != output & tmpout.exists()) {
					FileHelper.delete(tmpout, true);
				}
			}
			return null;
		}
	}
	/**
	 * Concatenates the input files in order.
	 * @param input input files.
	 * @param output output file
	 * @throws IOException
	 */
	public static void concat(SAMSequenceDictionary dictionary, List<File> input, File output) throws IOException {
		File tmpout = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(output, "gridss.tmp.concat.") : output;
		try (VariantContextWriter writer = new VariantContextWriterBuilder()
				.setOutputFile(tmpout)
				.setReferenceDictionary(dictionary)
				.unsetOption(Options.INDEX_ON_THE_FLY)
				.build()) {
			for (int i = 0; i < input.size(); i++) {
				try (VCFFileReader reader = new VCFFileReader(input.get(i), false)) {
					if (i == 0) {
						writer.writeHeader(reader.getFileHeader());
					}
					try (AsyncBufferedIterator<VariantContext> it = new AsyncBufferedIterator<>(reader.iterator(), input.get(i).getName())) {
						while (it.hasNext()) {
							writer.add(it.next());
						}
					}
				}
			}
			writer.close();
		}
		if (tmpout != output) {
			FileHelper.move(tmpout, output, true);
		}
	}
}
