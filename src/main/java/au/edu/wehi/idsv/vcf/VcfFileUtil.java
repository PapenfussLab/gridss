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
import java.util.Comparator;

import au.edu.wehi.idsv.ProcessingContext;

public class VcfFileUtil {
	private static final Log log = Log.getInstance(VcfFileUtil.class);
	/**
	 * Sorts a VCF according to the given sort order
	 * @param input unsorted input
	 * @param output sorted output to write
	 * @param sortComparator sort order
	 */
	public static void sort(ProcessingContext processContext, File input, File output, Comparator<VariantContext> sortComparator) {
		log.info("Sorting " + output);
		VCFFileReader reader = null;
		CloseableIterator<VariantContext> rit = null;
		VariantContextWriter writer = null;
		SortingCollection<VariantContext> collection = null;
		CloseableIterator<VariantContext> wit = null;
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
			collection.doneAdding();
			writer = processContext.getVariantContextWriter(output);
	    	wit = collection.iterator();
			while (wit.hasNext()) {
				writer.add(wit.next());
			}
		} finally {
			CloserUtil.close(writer);
			CloserUtil.close(wit);
			CloserUtil.close(rit);
			CloserUtil.close(reader);
			if (collection != null) collection.cleanup();
		}
	}
}
