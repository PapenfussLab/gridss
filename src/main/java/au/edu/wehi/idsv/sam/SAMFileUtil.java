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

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.concurrent.Callable;
import java.util.function.Function;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.IntermediateFileUtil;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.util.FileHelper;
import au.edu.wehi.idsv.validation.OrderAssertingIterator;

import com.google.common.collect.ImmutableList;

public class SAMFileUtil {
	private static final Log log = Log.getInstance(SAMFileUtil.class);
	/**
	 * Sorts records in the given SAM/BAM file by coordinate or queryname 
	 * @param unsorted input SAM/BAM file
	 * @param output sorted output file
	 * @param sortOrder sort order
	 * @throws IOException 
	 */
	public static void sort(ProcessingContext processContext, File unsorted, File output, SortOrder sortOrder) throws IOException {
		try {
			new SortCallable(processContext, unsorted, output, sortOrder, header -> header).call();
		} catch (IOException e) {
			log.error(log);
			throw new RuntimeException(e);
		}
	}
	/**
	 * Sorts records in the given SAM/BAM file by the given comparator 
	 * @param unsorted input SAM/BAM file
	 * @param output sorted output file
	 * @param sortComparator sort order
	 * @throws IOException 
	 */
	public static void sort(ProcessingContext processContext, File unsorted, File output, SAMRecordComparator sortComparator) {
		try {
		new SortCallable(processContext, unsorted, output, sortComparator, header -> header).call();
		} catch (IOException e) {
			log.error(log);
			throw new RuntimeException(e);
		}
	}
	/**
	 * SAM sorting task
	 * @author Daniel Cameron
	 *
	 */
	public static class SortCallable implements Callable<Void> {
		private final ProcessingContext processContext;
		private final File unsorted;
		private final File output;
		private final SAMRecordComparator sortComparator;
		private final SortOrder sortOrder;
		private Function<SAMFileHeader, SAMFileHeader> headerCallback;
		public SortCallable(ProcessingContext processContext, File unsorted, File output, SortOrder sortOrder, Function<SAMFileHeader, SAMFileHeader> headerCallback) {
			this(processContext, unsorted, output, sortOrder == SortOrder.coordinate ? new SAMRecordCoordinateComparator() : new SAMRecordQueryNameComparator(), sortOrder, headerCallback);
			switch (sortOrder) {
			case coordinate:
			case queryname:
				break;
			default:
				throw new IllegalArgumentException("Sort order not specified");
			}
		}
		public SortCallable(ProcessingContext processContext, File unsorted, File output, SAMRecordComparator sortComparator, Function<SAMFileHeader, SAMFileHeader> headerCallback) {
			this(processContext, unsorted, output, sortComparator, SortOrder.unsorted, headerCallback);
		}
		private SortCallable(ProcessingContext processContext, File unsorted, File output, SAMRecordComparator sortComparator, SortOrder sortOrder, Function<SAMFileHeader, SAMFileHeader> headerCallback) {
			this.processContext = processContext;
			this.unsorted = unsorted;
			this.output = output;
			this.sortComparator = sortComparator;
			this.sortOrder = sortOrder;
			this.headerCallback = headerCallback;
		}
		@Override
		public Void call() throws IOException {
			if (IntermediateFileUtil.checkIntermediate(output, unsorted, processContext.getConfig().ignoreFileTimestamps)) {
				log.info("Not sorting as output already exists: " + output);
				return null;
			}
			File tmpFile = FileSystemContext.getWorkingFileFor(output, "sorting");
			log.info("Sorting " + unsorted);
			SamReader reader = null;
			CloseableIterator<SAMRecord> rit = null;
			SAMFileWriter writer = null;
			CloseableIterator<SAMRecord> wit = null;
			SortingCollection<SAMRecord> collection = null;
			if (tmpFile.exists()) tmpFile.delete();
			try {
				reader = processContext.getSamReader(unsorted);
				SAMFileHeader header = reader.getFileHeader().clone();
				header.setSortOrder(sortOrder);
				if (headerCallback != null) {
					header = headerCallback.apply(header);
				}
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
				rit.close();
				rit = null;
				reader.close();
				reader = null;
				collection.doneAdding();
				writer = processContext.getSamFileWriterFactory(sortOrder == SortOrder.coordinate).makeSAMOrBAMWriter(header, true, tmpFile);
				writer.setProgressLogger(new ProgressLogger(log));
		    	wit = collection.iterator();
		    	if (Defaults.SANITY_CHECK_ITERATORS) {
					wit = new AutoClosingIterator<SAMRecord>(new OrderAssertingIterator<SAMRecord>(wit, sortComparator), ImmutableList.<Closeable>of(wit));
				}
				while (wit.hasNext()) {
					writer.addAlignment(wit.next());
				}
				wit.close();
				wit = null;
				writer.close();
				writer = null;
				collection.cleanup();
				collection = null;
				FileHelper.move(tmpFile, output, true);
			} finally {
				CloserUtil.close(writer);
				CloserUtil.close(wit);
				CloserUtil.close(rit);
				CloserUtil.close(reader);
				if (collection != null) collection.cleanup();
				if (tmpFile.exists()) tmpFile.delete();
			}
			return null;
		}
	}
}
