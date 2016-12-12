package au.edu.wehi.idsv.sam;

import java.io.File;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.concurrent.Callable;
import java.util.function.Function;

import com.google.common.collect.PeekingIterator;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.IntermediateFileUtil;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.FileHelper;
import au.edu.wehi.idsv.validation.OrderAssertingIterator;
import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordComparator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;

public class SAMFileUtil {
	private static final Log log = Log.getInstance(SAMFileUtil.class);
	/**
	 * Sorts records in the given SAM/BAM file by coordinate or queryname 
	 * @param unsorted input SAM/BAM file
	 * @param output sorted output file
	 * @param sortOrder sort order
	 * @throws IOException 
	 */
	public static void sort(FileSystemContext fsc, File unsorted, File output, SortOrder sortOrder) throws IOException {
		try {
			new SortCallable(fsc, unsorted, output, sortOrder, header -> header).call();
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
	public static void sort(FileSystemContext fsc, File unsorted, File output, SAMRecordComparator sortComparator) {
		try {
		new SortCallable(fsc, unsorted, output, sortComparator, header -> header).call();
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
		private final FileSystemContext fsc;
		private final File unsorted;
		private final File output;
		private final SAMRecordComparator sortComparator;
		private final SortOrder sortOrder;
		private final Function<SAMFileHeader, SAMFileHeader> headerCallback;
		private final SamReaderFactory readerFactory;
		private final SAMFileWriterFactory writerFactory;
		public SortCallable(FileSystemContext fsc, File unsorted, File output, SortOrder sortOrder, Function<SAMFileHeader, SAMFileHeader> headerCallback) {
			this(fsc, unsorted, output, null, sortOrder, headerCallback, null, null);
		}
		public SortCallable(FileSystemContext fsc, File unsorted, File output, SAMRecordComparator sortComparator, Function<SAMFileHeader, SAMFileHeader> headerCallback) {
			this(fsc, unsorted, output, sortComparator, SortOrder.unsorted, headerCallback, null, null);
		}
		private SortCallable(FileSystemContext fsc,
				File unsorted,
				File output,
				SAMRecordComparator sortComparator,
				SortOrder sortOrder,
				Function<SAMFileHeader, SAMFileHeader> headerCallback,
				SamReaderFactory readerFactory,
				SAMFileWriterFactory writerFactory) {
			this.fsc = fsc;
			this.unsorted = unsorted;
			this.output = output;
			this.sortComparator = sortComparator == null && sortOrder != null ? sortOrder.getComparatorInstance() : sortComparator;
			this.sortOrder = sortOrder;
			this.headerCallback = headerCallback;
			this.readerFactory = readerFactory == null ? SamReaderFactory.makeDefault() : readerFactory;
			this.writerFactory = writerFactory == null ? new SAMFileWriterFactory() : writerFactory;
			if (this.sortComparator == null) {
				throw new IllegalArgumentException("Sort order not specified");
			}
		}
		@Override
		public Void call() throws IOException {
			if (IntermediateFileUtil.checkIntermediate(output)) {
				log.info("Not sorting as output already exists: " + output);
				return null;
			}
			File tmpFile = FileSystemContext.getWorkingFileFor(output, "gridss.tmp.sorting.SAMFileUtil.");
			SortOrder existingSortOrder = getSortOrder(readerFactory, unsorted);  
			switch (existingSortOrder) {
				case coordinate:
				case queryname:
					if (sortOrder.equals(existingSortOrder)) {
						log.info(unsorted + " already sorted by " + sortOrder);
						FileHelper.copy(unsorted, tmpFile, true);
						FileHelper.move(tmpFile, output, true);
						return null;
					}
					break;
				default:
					break;
			}
			log.info("Sorting " + unsorted);
			SortingCollection<SAMRecord> collection = null;
			if (tmpFile.exists()) FileHelper.delete(tmpFile, true);
			try {
				SAMFileHeader header = null;
				try (SamReader reader = readerFactory.open(unsorted)) {
					header = reader.getFileHeader().clone();
					header.setSortOrder(sortOrder);
					if (headerCallback != null) {
						header = headerCallback.apply(header);
					}
					try (CloseableIterator<SAMRecord> rit = reader.iterator()) {
						collection = SortingCollection.newInstance(
								SAMRecord.class,
								new BAMRecordCodec(header),
								sortComparator,
								fsc.getMaxBufferedRecordsPerFile(),
								fsc.getTemporaryDirectory());
						while (rit.hasNext()) {
							collection.add(rit.next());
						}
					}
				}
				collection.doneAdding();
				try (SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(header, true, tmpFile)) {
					writer.setProgressLogger(new ProgressLogger(log, 10000000));
					try (CloseableIterator<SAMRecord> wit = collection.iterator()) {
						Iterator<SAMRecord> it = wit;
				    	if (Defaults.SANITY_CHECK_ITERATORS) {
							it = new OrderAssertingIterator<SAMRecord>(wit, sortComparator);
						}
						while (it.hasNext()) {
							writer.addAlignment(it.next());
						}
					}
				}
				collection.cleanup();
				collection = null;
				FileHelper.move(tmpFile, output, true);
			} finally {
				if (collection != null) collection.cleanup();
				if (tmpFile.exists()) FileHelper.delete(tmpFile, true);
			}
			return null;
		}
	}
	private static SortOrder getSortOrder(SamReaderFactory readerFactory, File file) throws IOException {
		try (SamReader reader = readerFactory.open(file)) {
			return reader.getFileHeader().getSortOrder();
		}
	}
	public static void merge(Collection<File> input, File output) throws IOException {
		merge(input, output, SamReaderFactory.makeDefault(), new SAMFileWriterFactory());
	}
	/**
	 * Merges a set of SAM files into a single file.
	 * The SAM header is taken from the first input file.
	 * @param input input files. All input files must have the same sort order.
	 * @param output output file.
	 * @param readerFactory
	 * @param writerFactory
	 * @throws IOException 
	 */
	public static void merge(Collection<File> input, File output, SamReaderFactory readerFactory, SAMFileWriterFactory writerFactory) throws IOException {
		if (input == null || input.size() == 0) return;
		File tmpFile = FileSystemContext.getWorkingFileFor(output, "gridss.tmp.merging.SAMFileUtil.");
		Map<SamReader, AsyncBufferedIterator<SAMRecord>> map = new HashMap<>(input.size());
		SAMFileHeader header = null;
		try {
			for (File in : input) {
				SamReader r = readerFactory.open(in);
				if (header == null) {
					header = r.getFileHeader();
				}
				map.put(r, new AsyncBufferedIterator<>(r.iterator(), in.getName()));
			}
			try (SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(header, true, tmpFile)) {
				Queue<PeekingIterator<SAMRecord>> queue = createMergeQueue(header.getSortOrder());
				for (PeekingIterator<SAMRecord> it : map.values()) {
					if (it.hasNext()) {
						queue.add(it);
					}
				}
				while (!queue.isEmpty()) {
					PeekingIterator<SAMRecord> it = queue.poll();
					SAMRecord r = it.next();
					writer.addAlignment(r);
					if (it.hasNext()) {
						queue.add(it);
					}
				}
			}
			for (Entry<SamReader, AsyncBufferedIterator<SAMRecord>> entry : map.entrySet()) {
				CloserUtil.close(entry.getValue());
				CloserUtil.close(entry.getKey());
			}
			FileHelper.move(tmpFile, output, true);
		} finally {
			for (Entry<SamReader, AsyncBufferedIterator<SAMRecord>> entry : map.entrySet()) {
				CloserUtil.close(entry.getValue());
				CloserUtil.close(entry.getKey());
			}
			if (tmpFile.exists()) {
				FileHelper.delete(tmpFile, true);
			}
		}
	}
	private static Queue<PeekingIterator<SAMRecord>> createMergeQueue(SortOrder sortOrder) {
		SAMRecordComparator comparator = sortOrder == null ? null : sortOrder.getComparatorInstance();
		if (comparator == null) return new ArrayDeque<>();
		return new PriorityQueue<>(Comparator.<PeekingIterator<SAMRecord>, SAMRecord>comparing(it -> it.peek(), comparator));
	}
}
