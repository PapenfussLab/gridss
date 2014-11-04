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
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.concurrent.Callable;

import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.IntermediateFileUtil;
import au.edu.wehi.idsv.ProcessingContext;

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
			new SortCallable(processContext, unsorted, output, sortOrder).call();
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
		new SortCallable(processContext, unsorted, output, sortComparator).call();
		} catch (IOException e) {
			log.error(log);
			throw new RuntimeException(e);
		}
	}
	/**
	 * SAM sorting task
	 * @author cameron.d
	 *
	 */
	public static class SortCallable implements Callable<Void> {
		private final ProcessingContext processContext;
		private final File unsorted;
		private final File output;
		private final SAMRecordComparator sortComparator;
		private final SortOrder sortOrder;
		public SortCallable(ProcessingContext processContext, File unsorted, File output, SortOrder sortOrder) {
			this(processContext, unsorted, output, sortOrder == SortOrder.coordinate ? new SAMRecordCoordinateComparator() : new SAMRecordQueryNameComparator(), SortOrder.unsorted);
			switch (sortOrder) {
			case coordinate:
			case queryname:
				break;
			default:
				throw new IllegalArgumentException("Sort order not specified");
			}
		}
		public SortCallable(ProcessingContext processContext, File unsorted, File output, SAMRecordComparator sortComparator) {
			this(processContext, unsorted, output, sortComparator, SortOrder.unsorted);
		}
		private SortCallable(ProcessingContext processContext, File unsorted, File output, SAMRecordComparator sortComparator, SortOrder sortOrder) {
			this.processContext = processContext;
			this.unsorted = unsorted;
			this.output = output;
			this.sortComparator = sortComparator;
			this.sortOrder = sortOrder;
		}
		@Override
		public Void call() throws IOException {
			if (IntermediateFileUtil.checkIntermediate(output, unsorted)) {
				log.info("Not sorting as output already exists: " + output);
				return null;
			}
			log.info("Sorting " + unsorted);
			SamReader reader = null;
			CloseableIterator<SAMRecord> rit = null;
			SAMFileWriter writer = null;
			CloseableIterator<SAMRecord> wit = null;
			SortingCollection<SAMRecord> collection = null;
			if (FileSystemContext.getWorkingFileFor(output).exists()) FileSystemContext.getWorkingFileFor(output).delete();
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
				writer = processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(header, true, FileSystemContext.getWorkingFileFor(output));
				writer.setProgressLogger(new ProgressLogger(log));
		    	wit = collection.iterator();
				while (wit.hasNext()) {
					writer.addAlignment(wit.next());
				}
				writer.close();
				collection.cleanup();
				Files.move(FileSystemContext.getWorkingFileFor(output).toPath(), output.toPath(), StandardCopyOption.REPLACE_EXISTING);
				File tmpIndex = new File(FileSystemContext.getWorkingFileFor(output).getAbsolutePath().substring(0, FileSystemContext.getWorkingFileFor(output).getAbsolutePath().length() - 4) + ".bai");
				if (tmpIndex.exists()) {
					File targetIndex = new File(output.getAbsolutePath().substring(0, output.getAbsolutePath().length() - 4) + ".bai");
					Files.move(tmpIndex.toPath(), targetIndex.toPath(), StandardCopyOption.REPLACE_EXISTING);
				}
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
