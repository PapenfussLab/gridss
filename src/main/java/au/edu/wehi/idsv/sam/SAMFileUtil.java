package au.edu.wehi.idsv.sam;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordComparator;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;

import java.io.File;

import au.edu.wehi.idsv.ProcessingContext;

public class SAMFileUtil {
	private static final Log log = Log.getInstance(SAMFileUtil.class);
	/**
	 * Sorts records in the given SAM/BAM file by coordinate or queryname 
	 * @param unsorted input SAM/BAM file
	 * @param output sorted output file
	 * @param sortOrder sort order
	 */
	public static void sort(ProcessingContext processContext, File unsorted, File output, SortOrder sortOrder) {
		switch (sortOrder) {
		case coordinate:
			sort(processContext, unsorted, output, new SAMRecordCoordinateComparator(), SortOrder.unsorted);
			break;
		case queryname:
			sort(processContext, unsorted, output, new SAMRecordQueryNameComparator(), SortOrder.unsorted);
			break;
		default:
			throw new IllegalArgumentException("Sort order not specified");
		}
	}
	/**
	 * Sorts records in the given SAM/BAM file by the given comparator 
	 * @param unsorted input SAM/BAM file
	 * @param output sorted output file
	 * @param sortComparator sort order
	 */
	public static void sort(ProcessingContext processContext, File unsorted, File output, SAMRecordComparator sortComparator) {
		sort(processContext, unsorted, output, sortComparator, SortOrder.unsorted);
	}
	private static void sort(ProcessingContext processContext, File unsorted, File output, SAMRecordComparator sortComparator, SortOrder sortOrder) {
		log.info("Sorting " + output);
		SamReader reader = null;
		CloseableIterator<SAMRecord> rit = null;
		SAMFileWriter writer = null;
		CloseableIterator<SAMRecord> wit = null;
		SortingCollection<SAMRecord> collection = null;
		try {
			reader = processContext.getSamReader(unsorted);
			SAMFileHeader header = reader.getFileHeader().clone();
			header.setSortOrder(sortOrder);
			rit = processContext.getSamReaderIterator(reader);
			collection = SortingCollection.newInstance(
					SAMRecord.class,
					new BAMRecordCodec(header),
					sortComparator,
					processContext.getFileSystemContext().getMaxBufferedRecordsPerFile(),
					processContext.getFileSystemContext().getTemporaryDirectory());
			while (rit.hasNext()) {
				collection.add(rit.next());
			}
			collection.doneAdding();
			writer = processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(header, true, output);
			writer.setProgressLogger(new ProgressLogger(log));
	    	wit = collection.iterator();
			while (wit.hasNext()) {
				writer.addAlignment(wit.next());
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
