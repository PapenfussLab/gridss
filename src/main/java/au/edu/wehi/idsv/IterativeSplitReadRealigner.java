package au.edu.wehi.idsv;

import au.edu.wehi.idsv.alignment.FastqAligner;
import au.edu.wehi.idsv.sam.NmTagIterator;
import au.edu.wehi.idsv.sam.SAMFileHeaderUtil;
import au.edu.wehi.idsv.sam.SAMFileUtil;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.FileHelper;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.PeekingIterator;
import com.google.common.io.Files;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.fastq.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class IterativeSplitReadRealigner extends SplitReadRealigner {
	private static final Log log = Log.getInstance(IterativeSplitReadRealigner.class);
	private final GenomicProcessingContext pc;
	private final FastqAligner aligner;
	private SamReaderFactory readerFactory;
	private SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
	private FastqWriterFactory fastqWriterFactory = new FastqWriterFactory();
	private List<File> tmpFiles = new ArrayList<>();

	public IterativeSplitReadRealigner(GenomicProcessingContext pc, FastqAligner aligner) {
		super(pc.getReference());
		this.pc = pc;
		this.readerFactory = SamReaderFactory.makeDefault().referenceSequence(pc.getReferenceFile());
		this.aligner = aligner;
	}

	public SamReaderFactory getReaderFactory() {
		return readerFactory;
	}
	public void setReaderFactory(SamReaderFactory readerFactory) {
		this.readerFactory = readerFactory;
	}
	public SAMFileWriterFactory getWriterFactory() {
		return writerFactory;
	}
	public void setWriterFactory(SAMFileWriterFactory writerFactory) {
		this.writerFactory = writerFactory;
	}
	public FastqWriterFactory getFastqWriterFactory() {
		return fastqWriterFactory;
	}
	public void setFastqWriterFactory(FastqWriterFactory fastqWriterFactory) {
		this.fastqWriterFactory = fastqWriterFactory;
	}

	@Override
	public void createSupplementaryAlignments(File input, File output, File unorderedOutput) throws IOException {
		try {
			int iteration = 0;
			File fq = pc.getFileSystemContext().getRealignmentFastq(input, iteration);
			File tmpfq = FileSystemContext.getWorkingFileFor(fq, "gridss.tmp.SplitReadRealigner.");
			int recordsWritten = createSupplementaryAlignmentFastq(input, tmpfq, false);
			Files.move(tmpfq, fq);
			tmpFiles.add(fq);
			tmpFiles.add(tmpfq);
			List<File> aligned = new ArrayList<>();
			while (recordsWritten > 0) {
				// Align
				File out = pc.getFileSystemContext().getRealignmentBam(input, iteration);
				File tmpout = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(out) : out;
				tmpFiles.add(out);
				tmpFiles.add(tmpout);
				aligner.align(fq, tmpout, pc.getReferenceFile(), getWorkerThreads(), getReference().getSequenceDictionary());
				if (tmpout != out) {
					FileHelper.move(tmpout, out, true);
				}
				aligned.add(out);
				// start next iteration
				iteration++;
				fq = pc.getFileSystemContext().getRealignmentFastq(out, iteration);
				tmpfq = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(fq) : fq;
				tmpFiles.add(fq);
				tmpFiles.add(tmpfq);
				recordsWritten = createSupplementaryAlignmentFastq(out, tmpfq, true);
				if (tmpfq != fq) {
					Files.move(tmpfq, fq);
				}
			}
			mergeSupplementaryAlignment(input, aligned, output, unorderedOutput);
		} finally {
			if (gridss.Defaults.DELETE_TEMPORARY_FILES) {
				for (File f : tmpFiles) {
					if (f.exists()) {
						FileHelper.delete(f, true);
					}
				}
			}
		}
	}
	public void mergeSupplementaryAlignment(File input, List<File> aligned, File output, File unorderedOutput) throws IOException {
		log.info("Merging split read alignments for ", output);
		File suppMerged = FileSystemContext.getWorkingFileFor(output, "gridss.tmp.SplitReadAligner.sa.");
		File tmpoutput = FileSystemContext.getWorkingFileFor(output);
		tmpFiles.add(suppMerged);
		tmpFiles.add(tmpoutput);
		List<SamReader> suppReaders = new ArrayList<>();
		List<PeekingIterator<SAMRecord>> suppIt = new ArrayList<>();
		SAMFileHeader header;
		try (SamReader reader = readerFactory.open(input)) {
			header = reader.getFileHeader();
			for (File sf : aligned) {
				SamReader suppReader = readerFactory.open(sf);
				suppReaders.add(suppReader);
				suppIt.add(new AsyncBufferedIterator<>(new NmTagIterator(suppReader.iterator(), getReference()), sf.getName()));
			}
			try (SAMFileWriter inputWriter = writerFactory.makeSAMOrBAMWriter(header, true, tmpoutput)) {
				SAMFileHeader suppUnsortedHeader = SAMFileHeaderUtil.minimal(header);
				suppUnsortedHeader.setSortOrder(SortOrder.unsorted);
				try (SAMFileWriter suppWriter = writerFactory.makeSAMOrBAMWriter(suppUnsortedHeader, true, suppMerged)) {
					try (AsyncBufferedIterator<SAMRecord> bufferedIt = new AsyncBufferedIterator<>(new NmTagIterator(reader.iterator(), getReference()), input.getName())) {
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
		if (unorderedOutput == null || output.equals(unorderedOutput)) {
			if (header.getSortOrder() != null && header.getSortOrder() != SortOrder.unsorted) {
				File suppMergedsorted = FileSystemContext.getWorkingFileFor(output, "gridss.tmp.SplitReadAligner.sorted.sa.");
				tmpFiles.add(suppMergedsorted);
				SAMFileUtil.sort(pc.getFileSystemContext(), suppMerged, suppMergedsorted, header.getSortOrder());
				FileHelper.move(suppMergedsorted, suppMerged, true);
			}
			SAMFileUtil.merge(ImmutableList.of(tmpoutput, suppMerged),
					output,
					pc.getSamReaderFactory(),
					pc.getSamFileWriterFactory());
		} else {
			FileHelper.move(tmpoutput, output, true);
			FileHelper.move(suppMerged, unorderedOutput, false);
		}
	}

	private void mergeSupplementaryAlignment(Iterator<SAMRecord> it, List<PeekingIterator<SAMRecord>> alignments, SAMFileWriter out, SAMFileWriter saout) {
		List<SAMRecord> salist = Lists.newArrayList();
		while (it.hasNext()) {
			salist.clear();
			SAMRecord r = it.next();
			String name = getEvidenceIdentifierGenerator().getAlignmentUniqueName(r);
			for (PeekingIterator<SAMRecord> sit : alignments) {
				while (sit.hasNext() && SplitReadHelper.getOriginatingAlignmentUniqueName(sit.peek()).equals(name)) {
					SAMRecord supp = sit.next();
					if (supp.getSupplementaryAlignmentFlag() || supp.isSecondaryAlignment()) {
						// only consider the best mapping location reported by the aligner
					} else {
						salist.add(supp);
					}
				}
			}
			writeCompletedAlignment(r, salist, out, saout);
		}
	}
	protected int createSupplementaryAlignmentFastq(File input, File fq, boolean isRecursive) throws IOException {
		int recordsWritten = 0;
		try (SamReader reader = readerFactory.open(input)) {
			try (AsyncBufferedIterator<SAMRecord> bufferedIt = new AsyncBufferedIterator<>(reader.iterator(), input.getName())) {
				try (FastqWriter writer = new AsyncFastqWriter(new BasicFastqWriter(fq), AsyncFastqWriter.DEFAULT_QUEUE_SIZE)) {
					while (bufferedIt.hasNext()) {
						SAMRecord r = bufferedIt.next();
						if (!shouldDropInputRecord(r)) {
							for (FastqRecord fqr : extract(r, isRecursive)) {
								writer.write(fqr);
								recordsWritten++;
							}
						}
					}
				}
			}
		}
		return recordsWritten;
	}
}
