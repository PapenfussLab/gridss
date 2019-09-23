package au.edu.wehi.idsv;

import au.edu.wehi.idsv.alignment.FastqAligner;
import au.edu.wehi.idsv.alignment.StreamingAligner;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.NmTagIterator;
import au.edu.wehi.idsv.sam.SAMFileUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
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
import java.util.*;

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
	private boolean adjustPrimary = false;
	private ReferenceLookup reference = null;
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
	public void createSupplementaryAlignments(final StreamingAligner aligner, final File input, final File output, final File outputModified, final boolean rewriteOA, final int maxBufferedRecords) throws IOException {
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
		int recordNumber = 0;
		try (SamReader reader = readerFactory.open(input)) {
			SAMFileHeader header = reader.getFileHeader().clone();
			header.setSortOrder(SortOrder.unsorted);
			SAMFileWriter modifiedWriter = null;
			try (SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(header, true, output)) {
				modifiedWriter =  writer;
				if (outputModified != null && !output.equals(outputModified)) {
					modifiedWriter = writerFactory.makeSAMOrBAMWriter(header, true, outputModified);
				}
				try (AsyncBufferedIterator<SAMRecord> bufferedIt = new AsyncBufferedIterator<>(reader.iterator(), input.getName())) {
					while (bufferedIt.hasNext()) {
						SAMRecord r = bufferedIt.next();
						processInputRecord(aligner, rootExtractor, realignments, writer, r);
						if (aligner.outstandingAlignmentRecord() >= maxBufferedRecords) {
							log.info(String.format("%d records awaiting alignment by external aligner. Flushing.", maxBufferedRecords));
							aligner.flush();
						}
						while (aligner.hasAlignmentRecord()) {
							processAlignmentRecord(aligner, recursiveExtractor, realignments, writer, modifiedWriter, rewriteOA);
						}
						if (++recordNumber % 1000 == 0) {
							String msg = String.format("Processed %d records. %d in aligner input buffer. %d in aligner output buffer. %s records in lookup", recordNumber, aligner.outstandingAlignmentRecord(), aligner.processedAlignmentRecords(), realignments.size());
							log.debug(msg);
							if (recordNumber % 1000000 == 0) {
								log.info(msg);
							}
						}
					}
					// flush out all realignments
					aligner.flush();
					while (aligner.hasAlignmentRecord()) {
						// perform nested realignment
						while (aligner.hasAlignmentRecord()) {
							processAlignmentRecord(aligner, recursiveExtractor, realignments, writer, modifiedWriter, rewriteOA);
						}
						aligner.flush();
					}
				}
				if (realignments.size() != 0) {
					log.error(String.format("External aligner did not return alignments for %d records including %s.", realignments.size(), realignments.values().iterator().next().originatingRecord.getReadName()));
					for (SplitReadRealignmentInfo info : realignments.values()) {
						writeCompletedAlignment(info.originatingRecord, info.realignments, writer, modifiedWriter, rewriteOA);
					}
				}
			} finally {
				CloserUtil.close(modifiedWriter);
			}
		}
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
			SAMFileWriter writer, SAMFileWriter modifiedRecordWriter, boolean writeOA) throws IOException {
		SAMRecord supp = aligner.getAlignment();
		String lookupkey = SplitReadHelper.getOriginatingAlignmentUniqueName(supp);
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
				writeCompletedAlignment(info.originatingRecord, info.realignments, writer, modifiedRecordWriter, writeOA);
				realignments.remove(lookupkey);
			} else {
				//log.debug(String.format("%s: %d outstanding alignments", info.originatingRecord.getReadName(), info.outstandingRealignments));
			}
		}
	}
	private void writeCompletedAlignment(SAMRecord primary, List<SAMRecord> realignments, SAMFileWriter recordWriter, SAMFileWriter realignmentWriter, boolean writeOA) {
		if (isRealignExistingSplitReads() || isRealignEntireRecord()) {
			if (primary.getSupplementaryAlignmentFlag()) {
				// If we're realigning, we need to drop all existing supplementary alignments
				return;
			}
			primary.setAttribute(SAMTag.SA.name(), null);
		}
		boolean primaryHasMoved = prepareRecordsForWriting(primary, realignments, writeOA);
		if (primary.getReadUnmappedFlag() || primaryHasMoved) {
			// we'll break sort ordering if we write it back to the input file
			realignmentWriter.addAlignment(primary);
		} else {
			recordWriter.addAlignment(primary);
		}
		for (SAMRecord sar : realignments) {
			realignmentWriter.addAlignment(sar);
		}
	}
	private boolean prepareRecordsForWriting(SAMRecord primary, List<SAMRecord> realignments, boolean writeOA) {
		int primaryReferenceIndex = primary.getReadUnmappedFlag() ? -1 : primary.getReferenceIndex();
		int primaryAlignmentStart = primary.getReadUnmappedFlag() ? -1 : primary.getAlignmentStart();
		// Only do anything to records that we actually attempted realignment for
		if (realignments.size() > 0) {
			if (isRealignEntireRecord() &&
					// special case exclusion of unanchored assemblies as
					// we repurpose the CIGAR string to encode the breakend interval
					// and realigning will break that
					!AssemblyAttributes.isUnanchored(primary)) {
				SAMRecord newPrimaryAlignmentPosition = SplitReadHelper.replaceAlignment(primary, realignments, writeOA);
				if (newPrimaryAlignmentPosition != null && !realignments.remove(newPrimaryAlignmentPosition)) {
					throw new RuntimeException("Sanity check failure: no supplementary alignment was removed when replacing alignment");
				}
			}
			SplitReadHelper.convertToSplitRead(primary, realignments, reference, isAdjustPrimaryAlignment() || isRealignEntireRecord());
			for (int i = 0; i < realignments.size(); i++) {
				SAMRecord r = realignments.get(i);
				if (r.getReadUnmappedFlag() || SAMRecordUtil.getStartClipLength(r) == SAMRecordUtil.getReadLengthIncludingHardClipping(r)) {
					realignments.remove(i);
					i--;
				}
			}
			// TODO: remove alignments which are contained by another alignment
			if (primary.getReadUnmappedFlag() || SAMRecordUtil.getStartClipLength(primary) == SAMRecordUtil.getReadLengthIncludingHardClipping(primary)) {
				SplitReadHelper.replaceAlignment(primary, realignments, writeOA);
			}
		}
		return primary.getReadUnmappedFlag() || primary.getReferenceIndex() != primaryReferenceIndex || primary.getAlignmentStart() != primaryAlignmentStart;
	}
	public void createSupplementaryAlignments(FastqAligner aligner, File input, File output, File unorderedOutput, final boolean rewriteOA) throws IOException {
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
			mergeSupplementaryAlignment(input, aligned, output, unorderedOutput, rewriteOA);
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
	private void mergeSupplementaryAlignment(File input, List<File> aligned, File output, File unorderedOutput, boolean rewriteOA) throws IOException {
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
						mergeSupplementaryAlignment(bufferedIt, suppIt, inputWriter, suppWriter, rewriteOA);
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
					pc.getSamFileWriterFactory(header.getSortOrder() == SortOrder.coordinate));
		} else {
			FileHelper.move(tmpoutput, output, true);
			FileHelper.move(suppMerged, unorderedOutput, false);
		}
	}
	private void mergeSupplementaryAlignment(Iterator<SAMRecord> it, List<PeekingIterator<SAMRecord>> alignments, SAMFileWriter out, SAMFileWriter saout, boolean writeOA) {
		List<SAMRecord> salist = Lists.newArrayList();
		while (it.hasNext()) {
			salist.clear();
			SAMRecord r = it.next();
			String name = eidgen.getAlignmentUniqueName(r);
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
			writeCompletedAlignment(r, salist, out, saout, writeOA);
		}
	}
	protected int createSupplementaryAlignmentFastq(File input, File fq, boolean isRecursive) throws IOException {
		int recordsWritten = 0;
		try (SamReader reader = readerFactory.open(input)) {
			try (AsyncBufferedIterator<SAMRecord> bufferedIt = new AsyncBufferedIterator<>(reader.iterator(), input.getName())) {
				try (FastqWriter writer = new AsyncFastqWriter(new NonFlushingBasicFastqWriter(fq), AsyncFastqWriter.DEFAULT_QUEUE_SIZE)) {
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
	public ReferenceLookup getReference() {
		return reference;
	}
	public void setReference(ReferenceLookup reference) {
		this.reference = reference;
	}
	public void setAdjustPrimaryAlignment(boolean adjustPrimary) {
		this.adjustPrimary = adjustPrimary;
	}
	public boolean isAdjustPrimaryAlignment() {
		return adjustPrimary;
	}
}
