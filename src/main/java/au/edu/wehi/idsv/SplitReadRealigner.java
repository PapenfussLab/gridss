package au.edu.wehi.idsv;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.PeekingIterator;
import com.google.common.io.Files;

import au.edu.wehi.idsv.alignment.FastqAligner;
import au.edu.wehi.idsv.sam.NmTagIterator;
import au.edu.wehi.idsv.sam.SAMFileUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.util.CloserUtil;

public class SplitReadRealigner {
	private GenomicProcessingContext pc;
	private int minSoftClipLength = 15;
	private SamReaderFactory readerFactory = SamReaderFactory.make();
	private SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
	private FastqWriterFactory fastqWriterFactory = new FastqWriterFactory();
	private boolean processSecondaryAlignments = false;
	private int workerThreads = Runtime.getRuntime().availableProcessors();
	private FastqAligner aligner;
	
	public SplitReadRealigner(GenomicProcessingContext pc, FastqAligner aligner) {
		this.pc = pc;
		this.aligner = aligner;
	}
	public void createSupplementaryAlignments(File input, File output) throws IOException {
		int iteration = 0;
		File fq = pc.getFileSystemContext().getRealignmentFastq(input, iteration);
		File tmpfq = FileSystemContext.getWorkingFileFor(fq);
		int recordsWritten = createSupplementaryAlignmentFastq(input, tmpfq, false);
		Files.move(tmpfq, fq);
		List<File> aligned = new ArrayList<>();
		while (recordsWritten > 0) {
			// Align
			File out = pc.getFileSystemContext().getRealignmentBam(input, iteration);
			File tmpout = FileSystemContext.getWorkingFileFor(out);
			aligner.align(fq, tmpout, pc.getReferenceFile(), workerThreads);
			Files.move(tmpout, out);
			aligned.add(out);
			// start next iteration
			iteration++;
			fq = pc.getFileSystemContext().getRealignmentFastq(out, iteration);
			tmpfq = FileSystemContext.getWorkingFileFor(fq);
			recordsWritten = createSupplementaryAlignmentFastq(out, tmpfq, true);
			Files.move(tmpfq, fq);
		}
		mergeSupplementaryAlignment(input, aligned, output);
	}
	private void mergeSupplementaryAlignment(File input, List<File> aligned, File output) throws IOException {
		File suppMerged = FileSystemContext.getWorkingFileFor(output, "tmp.sa.");
		File tmpoutput = FileSystemContext.getWorkingFileFor(output);
		List<SamReader> suppReaders = new ArrayList<>();
		List<PeekingIterator<SAMRecord>> suppIt = new ArrayList<>();
		SAMFileHeader header;
		try (SamReader reader = readerFactory.open(input)) {
			header = reader.getFileHeader();
			for (File sf : aligned) {
				SamReader suppReader = readerFactory.open(sf);
				suppReaders.add(suppReader);
				suppIt.add(new AsyncBufferedIterator<>(new NmTagIterator(suppReader.iterator(), pc.getReference()), sf.getName()));
			}
			try (SAMFileWriter inputWriter = writerFactory.makeSAMOrBAMWriter(header, true, tmpoutput)) {
				try (SAMFileWriter suppWriter = writerFactory.makeSAMOrBAMWriter(header, false, suppMerged)) {
					try (AsyncBufferedIterator<SAMRecord> bufferedIt = new AsyncBufferedIterator<>(new NmTagIterator(reader.iterator(), pc.getReference()), input.getName())) {
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
		if (header.getSortOrder() != null && header.getSortOrder() != SortOrder.unsorted) {
			File suppMergedsorted = FileSystemContext.getWorkingFileFor(output, "tmp.sorted.sa.");
			SAMFileUtil.sort(pc.getFileSystemContext(), suppMerged, suppMergedsorted, header.getSortOrder());
			Files.move(suppMergedsorted, suppMerged);
		}
		SAMFileUtil.merge(ImmutableList.of(tmpoutput, suppMerged), output);
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
			try (AsyncBufferedIterator<SAMRecord> bufferedIt = new AsyncBufferedIterator<>(reader.iterator(), input.getName())) {
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
