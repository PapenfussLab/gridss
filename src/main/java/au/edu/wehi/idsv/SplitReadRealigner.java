package au.edu.wehi.idsv;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.lang.NotImplementedException;

import com.google.common.collect.Lists;
import com.google.common.collect.PeekingIterator;
import com.google.common.io.Files;

import au.edu.wehi.idsv.alignment.FastqAligner;
import au.edu.wehi.idsv.picard.SynchronousReferenceLookupAdapter;
import au.edu.wehi.idsv.sam.NmTagIterator;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.CloserUtil;

public class SplitReadRealigner {
	private static final int ASYNC_BUFFERS = 2;
	private static final int ASYNC_BUFFER_SIZE = 300;
	private FileSystemContext fsc;
	private int minSoftClipLength = 15;
	private SamReaderFactory readerFactory = SamReaderFactory.make();
	private SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
	private FastqWriterFactory fastqWriterFactory = new FastqWriterFactory();
	private boolean processSecondaryAlignments = false;
	private int workerThreads = Runtime.getRuntime().availableProcessors();
	private FastqAligner aligner;
	
	public SplitReadRealigner(FileSystemContext fsc, FastqAligner aligner) {
		this.fsc = fsc;
		this.aligner = aligner;
	}
	public void createSupplementaryAlignments(File input, File referenceFile, File output) throws IOException {
		createSupplementaryAlignments(input, referenceFile, new SynchronousReferenceLookupAdapter(new IndexedFastaSequenceFile(referenceFile)), output);
	}
	public void createSupplementaryAlignments(File input, File referenceFile, SynchronousReferenceLookupAdapter reference, File output) throws IOException {
		int iteration = 0;
		File fq = fsc.getRealignmentFastq(input, iteration);
		File tmpfq = FileSystemContext.getWorkingFileFor(fq);
		int recordsWritten = createSupplementaryAlignmentFastq(input, tmpfq, false);
		Files.move(tmpfq, fq);
		List<File> aligned = new ArrayList<>();
		while (recordsWritten > 0) {
			// Align
			File out = fsc.getRealignmentBam(input, iteration);
			File tmpout = FileSystemContext.getWorkingFileFor(out);
			aligner.align(fq, tmpout, referenceFile, workerThreads);
			Files.move(tmpout, out);
			aligned.add(out);
			// start next iteration
			iteration++;
			fq = fsc.getRealignmentFastq(input, iteration);
			tmpfq = FileSystemContext.getWorkingFileFor(fq);
			recordsWritten = createSupplementaryAlignmentFastq(out, tmpfq, true);
			Files.move(tmpfq, fq);
		}
		mergeSupplementaryAlignment(input, aligned, output, reference);
	}
	private void mergeSupplementaryAlignment(File input, List<File> aligned, File output, ReferenceSequenceFile reference) throws IOException {
		File suppMerged = FileSystemContext.getWorkingFileFor(output, "tmp.sa.");
		File tmpoutput = FileSystemContext.getWorkingFileFor(output);
		List<SamReader> suppReaders = new ArrayList<>();
		List<PeekingIterator<SAMRecord>> suppIt = new ArrayList<>();
		try (SamReader reader = readerFactory.open(input)) {
			SAMFileHeader header = reader.getFileHeader();
			for (File sf : aligned) {
				SamReader suppReader = readerFactory.open(sf);
				suppReaders.add(suppReader);
				suppIt.add(new AsyncBufferedIterator<>(new NmTagIterator(suppReader.iterator(), reference), ASYNC_BUFFERS, ASYNC_BUFFER_SIZE));
			}
			try (SAMFileWriter inputWriter = writerFactory.makeSAMOrBAMWriter(header, true, tmpoutput)) {
				try (SAMFileWriter suppWriter = writerFactory.makeSAMOrBAMWriter(header, false, suppMerged)) {
					try (AsyncBufferedIterator<SAMRecord> bufferedIt = new AsyncBufferedIterator<>(reader.iterator(), ASYNC_BUFFERS, ASYNC_BUFFER_SIZE)) {
						mergeSupplementaryAlignment(bufferedIt, suppIt, inputWriter, suppWriter);
					}
				}
			}
		} finally {
			for (Iterator<SAMRecord> it : suppIt) {
				CloserUtil.close(it);
			}
			for (SamReader sr : suppReaders) {
				sr.close();
			}
		}
		throw new NotImplementedException();
		// sort suppMerged (by input header order)
		// merge suppMerged tmpoutput > output
	}
	private void mergeSupplementaryAlignment(Iterator<SAMRecord> it, List<PeekingIterator<SAMRecord>> alignments, SAMFileWriter out, SAMFileWriter saout) {
		List<SAMRecord> salist = Lists.newArrayList();
		while (it.hasNext()) {
			salist.clear();
			SAMRecord r = it.next();
			String name = SAMRecordUtil.getAlignmentUniqueName(r);
			for (PeekingIterator<SAMRecord> sit : alignments) {
				while (sit.hasNext() && SplitReadIdentificationHelper.getAlignmentUniqueName(sit.peek()).equals(name)) {
					salist.add(sit.next());
				}
			}
			if (salist.size() > 0) {
				SplitReadIdentificationHelper.convertToSplitRead(r, salist);
			}
			out.addAlignment(r);
			for (SAMRecord sar : salist) {
				if (!sar.getReadUnmappedFlag()) {
					saout.addAlignment(sar);
				}
			}
		}
	}
	private int createSupplementaryAlignmentFastq(File input, File fq, boolean isRecursive) throws IOException {
		int recordsWritten = 0;
		try (SamReader reader = readerFactory.open(input)) {
			try (AsyncBufferedIterator<SAMRecord> bufferedIt = new AsyncBufferedIterator<>(reader.iterator(), ASYNC_BUFFERS, ASYNC_BUFFER_SIZE)) {
				try (FastqWriter writer = fastqWriterFactory.newWriter(fq)) {
					SplitReadFastqExtractionIterator fastqit = new SplitReadFastqExtractionIterator(bufferedIt, isRecursive, minSoftClipLength, processSecondaryAlignments);
					while (fastqit.hasNext()) {
						writer.write(fastqit.next());
						recordsWritten++;
					}
					
				}
			}
		}
		return recordsWritten;
	}
}
