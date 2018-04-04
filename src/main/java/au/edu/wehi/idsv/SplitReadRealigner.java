package au.edu.wehi.idsv;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.PeekingIterator;
import com.google.common.io.Files;

import au.edu.wehi.idsv.alignment.FastqAligner;
import au.edu.wehi.idsv.alignment.StreamingAligner;
import au.edu.wehi.idsv.sam.NmTagIterator;
import au.edu.wehi.idsv.sam.SAMFileUtil;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.FileHelper;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

public class SplitReadRealigner {
	private static final Log log = Log.getInstance(SplitReadRealigner.class);
	private final GenomicProcessingContext pc;
	private int minSoftClipLength = 1;
	private float minSoftClipQuality = 0;
	private int workerThreads = Runtime.getRuntime().availableProcessors();
	private SamReaderFactory readerFactory;
	private SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
	private FastqWriterFactory fastqWriterFactory = new FastqWriterFactory();
	private boolean processSecondaryAlignments = false;
	private boolean realignExistingSplitReads = false;
	private boolean realignEntireRecord = false;
	/**
	 * Alignment-unique read identifier must be hashed to ensure that the read names
	 * written to the fastq files do not exceed the BAM limit of 254 even when the
	 * input reads names are at this limit.
	 * See https://github.com/PapenfussLab/gridss/issues/82
	 */
	private EvidenceIdentifierGenerator eidgen = new HashedEvidenceIdentifierGenerator(); 
	private List<File> tmpFiles = new ArrayList<>();
	
	public SplitReadRealigner(GenomicProcessingContext pc) {
		this.pc = pc;
		this.readerFactory = SamReaderFactory.makeDefault().referenceSequence(pc.getReferenceFile());
	}
	public int getMinSoftClipLength() {
		return minSoftClipLength;
	}
	public void setMinSoftClipLength(int minSoftClipLength) {
		this.minSoftClipLength = minSoftClipLength;
	}
	public float getMinSoftClipQuality() {
		return minSoftClipQuality;
	}
	public void setMinSoftClipQuality(float minSoftClipQuality) {
		this.minSoftClipQuality = minSoftClipQuality;
	}
	public int getWorkerThreads() {
		return workerThreads;
	}
	public void setWorkerThreads(int workerThreads) {
		this.workerThreads = workerThreads;
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
	public boolean isProcessSecondaryAlignments() {
		return processSecondaryAlignments;
	}
	public void setProcessSecondaryAlignments(boolean processSecondaryAlignments) {
		this.processSecondaryAlignments = processSecondaryAlignments;
	}
	private class SplitReadRealignmentInfo {
		@Override
		public int hashCode() {
			return lookupKey.hashCode();
		}
		@Override
		public boolean equals(Object obj) {
			return lookupKey.equals(((SplitReadRealignmentInfo)obj).lookupKey);
		}
		public SplitReadRealignmentInfo(SAMRecord record) {
			this.originatingRecord = record;
			this.lookupKey = eidgen.getAlignmentUniqueName(record);
		}
		public SAMRecord originatingRecord;
		public String lookupKey;
		public List<SAMRecord> realignments = new ArrayList<>(2);
		public int outstandingRealignments = 0;
	}
	public void createSupplementaryAlignments(StreamingAligner aligner, File input, File output) throws IOException {
		SplitReadFastqExtractor rootExtractor = new SplitReadFastqExtractor(false,
				minSoftClipLength,
				minSoftClipQuality,
				isProcessSecondaryAlignments(),
				isRealignExistingSplitReads(),
				isRealignEntireRecord(),
				eidgen);
		SplitReadFastqExtractor recursiveExtractor = new SplitReadFastqExtractor(true,
				minSoftClipLength,
				minSoftClipQuality,
				false,
				isRealignExistingSplitReads(),
				isRealignEntireRecord(),
				eidgen);
		
		Map<String, SplitReadRealignmentInfo> realignments = new HashMap<>();
		
		try (SamReader reader = readerFactory.open(input)) {
			SAMFileHeader header = reader.getFileHeader().clone();
			header.setSortOrder(SortOrder.unsorted);
			try (SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(header, true, output)) {
				try (AsyncBufferedIterator<SAMRecord> bufferedIt = new AsyncBufferedIterator<>(reader.iterator(), input.getName())) {
					while (bufferedIt.hasNext()) {
						SAMRecord r = bufferedIt.next();
						processInputRecord(aligner, rootExtractor, realignments, writer, r);
						while (aligner.hasAlignmentRecord()) {
							processAlignmentRecord(aligner, recursiveExtractor, realignments, writer);
						}
					}
					// flush out all realignments
					aligner.flush();
					while (aligner.hasAlignmentRecord()) {
						// perform nested realignment
						while (aligner.hasAlignmentRecord()) {
							processAlignmentRecord(aligner, recursiveExtractor, realignments, writer);
						}
						aligner.flush();
					}
				}
			}
		}
		assert(realignments.size() == 0);
	}
	private void processInputRecord(StreamingAligner aligner, SplitReadFastqExtractor rootExtractor,
			Map<String, SplitReadRealignmentInfo> realignments, SAMFileWriter writer, SAMRecord r) throws IOException {
		if ((realignEntireRecord || realignExistingSplitReads) && r.getSupplementaryAlignmentFlag()) {
			// drop existing supp alignments
			return;
		}
		List<FastqRecord> softclipRealignments = rootExtractor.extract(r);
		if (softclipRealignments.size() == 0) {
			// nothing to do - just output the record
			writer.addAlignment(r);
		} else {
			// perform split read realignment
			SplitReadRealignmentInfo info = new SplitReadRealignmentInfo(r);
			realignments.put(info.lookupKey, info);
			for (FastqRecord fq : softclipRealignments) {
				aligner.asyncAlign(fq);
				info.outstandingRealignments++;
			}
		}
	}
	private void processAlignmentRecord(StreamingAligner aligner,
			SplitReadFastqExtractor recursiveExtractor, Map<String, SplitReadRealignmentInfo> realignments,
			SAMFileWriter writer) throws IOException {
		SAMRecord supp = aligner.getAlignment();
		String lookupkey = SplitReadIdentificationHelper.getOriginatingAlignmentUniqueName(supp);
		SplitReadRealignmentInfo info = realignments.get(lookupkey);
		if (supp.getSupplementaryAlignmentFlag() || supp.isSecondaryAlignment()) {
			// only consider the best mapping location reported by the aligner
			//log.debug(String.format("%s: ignoring supp alignment", supp.getReadName()));
		} else {
			assert(info.outstandingRealignments > 0);
			info.outstandingRealignments--;
			if (!supp.getReadUnmappedFlag()) {
				info.realignments.add(supp);
				List<FastqRecord> nestedRealignments = recursiveExtractor.extract(supp);
				for (FastqRecord fq : nestedRealignments) {
					aligner.asyncAlign(fq);
					info.outstandingRealignments++;
					//log.debug(String.format("%s: performing nested realignment. %d realignments now outstanding", info.originatingRecord.getReadName(), info.outstandingRealignments));
				}
			}
			// all splits identified
			if (info.outstandingRealignments == 0) {
				writeCompletedAlignment(info.originatingRecord, info.realignments, writer, writer);
				realignments.remove(lookupkey);
			} else {
				//log.debug(String.format("%s: %d outstanding alignments", info.originatingRecord.getReadName(), info.outstandingRealignments));
			}
		}
	}
	private void writeCompletedAlignment(SAMRecord record, List<SAMRecord> realignments, SAMFileWriter recordWriter, SAMFileWriter realignmentWriter) {
		if (isRealignExistingSplitReads() || isRealignEntireRecord()) {
			if (record.getSupplementaryAlignmentFlag()) {
				// If we're realigning, we need to drop all existing supplementary alignments
				return;
			}
			record.setAttribute(SAMTag.SA.name(), null);
		}
		// Only do anything to records that we actually attempted realignment for 
		if (realignments.size() > 0) {
			if (isRealignEntireRecord() &&
					// special case exclusion of unanchored assemblies as
					// we repurpose the CIGAR string to encode the breakend interval
					// and realigning will break that
					!AssemblyAttributes.isUnanchored(record)) {
				SAMRecord newPrimaryAlignmentPosition = SplitReadIdentificationHelper.replaceAlignment(record, realignments);
				if (newPrimaryAlignmentPosition != null && !realignments.remove(newPrimaryAlignmentPosition)) {
					throw new RuntimeException("Sanity check failure: no supplementary alignment was removed when replacing alignment");
				}
			}
			SplitReadIdentificationHelper.convertToSplitRead(record, realignments);
		}
		recordWriter.addAlignment(record);
		for (SAMRecord sar : realignments) {
			realignmentWriter.addAlignment(sar);
		}
	}
	public void createSupplementaryAlignments(FastqAligner aligner, File input, File output) throws IOException {
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
				aligner.align(fq, tmpout, pc.getReferenceFile(), workerThreads);
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
			mergeSupplementaryAlignment(input, aligned, output);
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
	private void mergeSupplementaryAlignment(File input, List<File> aligned, File output) throws IOException {
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
			File suppMergedsorted = FileSystemContext.getWorkingFileFor(output, "gridss.tmp.SplitReadAligner.sorted.sa.");
			tmpFiles.add(suppMergedsorted);
			SAMFileUtil.sort(pc.getFileSystemContext(), suppMerged, suppMergedsorted, header.getSortOrder());
			FileHelper.move(suppMergedsorted, suppMerged, true);
		}
		SAMFileUtil.merge(ImmutableList.of(tmpoutput, suppMerged),
				output,
				pc.getSamReaderFactory(),
				pc.getSamFileWriterFactory(header.getSortOrder() == SortOrder.coordinate));
	}
	private void mergeSupplementaryAlignment(Iterator<SAMRecord> it, List<PeekingIterator<SAMRecord>> alignments, SAMFileWriter out, SAMFileWriter saout) {
		List<SAMRecord> salist = Lists.newArrayList();
		while (it.hasNext()) {
			salist.clear();
			SAMRecord r = it.next();
			String name = eidgen.getAlignmentUniqueName(r);
			for (PeekingIterator<SAMRecord> sit : alignments) {
				while (sit.hasNext() && SplitReadIdentificationHelper.getOriginatingAlignmentUniqueName(sit.peek()).equals(name)) {
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
				try (FastqWriter writer = fastqWriterFactory.newWriter(fq)) {
					SplitReadFastqExtractionIterator fastqit = new SplitReadFastqExtractionIterator(
							bufferedIt,
							isRecursive,
							minSoftClipLength,
							minSoftClipQuality,
							!isRecursive && isProcessSecondaryAlignments(),
							isRealignExistingSplitReads(),
							isRealignEntireRecord(),
							eidgen);
					while (fastqit.hasNext()) {
						writer.write(fastqit.next());
						recordsWritten++;
					}
					
				}
			}
		}
		return recordsWritten;
	}
	public boolean isRealignExistingSplitReads() {
		return realignExistingSplitReads;
	}
	public void setRealignExistingSplitReads(boolean realignExistingSplitReads) {
		this.realignExistingSplitReads = realignExistingSplitReads;
	}
	public boolean isRealignEntireRecord() {
		return realignEntireRecord;
	}
	public void setRealignEntireRecord(boolean realignEntireRecord) {
		this.realignEntireRecord = realignEntireRecord;
	}
}
