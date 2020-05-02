package au.edu.wehi.idsv;

import au.edu.wehi.idsv.alignment.StreamingAligner;
import au.edu.wehi.idsv.sam.SAMFileHeaderUtil;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import htsjdk.samtools.*;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class StreamingSplitReadRealigner extends SplitReadRealigner {
    private static final Log log = Log.getInstance(StreamingSplitReadRealigner.class);
    private final StreamingAligner aligner;
    private final GenomicProcessingContext pc;
    private final int maxBufferedRecords;

    public StreamingSplitReadRealigner(GenomicProcessingContext pc, StreamingAligner aligner, int maxBufferedRecords) {
        super(pc.getReference());
        this.pc = pc;
        this.aligner = aligner;
        this.maxBufferedRecords = maxBufferedRecords;
    }

    public void process(Iterator<SAMRecord> it, SAMFileWriter coordinateSortedWriter, SAMFileWriter unorderedWriter) throws IOException {
        Map<String, SplitReadRealignmentInfo> lookup = new HashMap<>();
        ProgressLogger progress = new ProgressLogger(log);
        int recordNumber = 0;
        while (it.hasNext()) {
            if (++recordNumber % 1000 == 0) {
                String msg = String.format("Processed %d records. %d in aligner input buffer. %d in aligner output buffer. %s records in lookup", recordNumber, aligner.outstandingAlignmentRecord(), aligner.processedAlignmentRecords(), lookup.size());
                log.debug(msg);
                if (recordNumber % 1000000 == 0) {
                    log.info(msg);
                }
            }
            SAMRecord r = it.next();
            progress.record(r);
            processCompletedAsyncRealignments(lookup, coordinateSortedWriter, unorderedWriter);
            processInputRecord(r, lookup, coordinateSortedWriter);
        }
        // perform nested realignment to ensure all records are fully recursively realigned
        aligner.flush();
        while (aligner.processedAlignmentRecords() > 0) {
            processCompletedAsyncRealignments(lookup, coordinateSortedWriter, unorderedWriter);
            aligner.flush();
        }
    }

    private void processCompletedAsyncRealignments(
            Map<String, SplitReadRealignmentInfo> lookup,
            SAMFileWriter coordinateSortedWriter,
            SAMFileWriter unorderedWriter) throws IOException {
        if (aligner.processedAlignmentRecords() > 0) {
            SAMRecord realignment = aligner.getAlignment();
            processAlignmentRecord(realignment, lookup, coordinateSortedWriter, unorderedWriter);
        }
        if (aligner.outstandingAlignmentRecord() >= maxBufferedRecords) {
            log.info(String.format("%d records awaiting alignment by external aligner. Flushing.", maxBufferedRecords));
            aligner.flush();
            processCompletedAsyncRealignments(lookup, coordinateSortedWriter, unorderedWriter);
        }
    }

    private void processInputRecord(SAMRecord record, Map<String, SplitReadRealignmentInfo> realignments, SAMFileWriter coordinateSortedWriter) throws IOException {
        if (shouldDropInputRecord(record)) {
            return;
        }
        List<FastqRecord> softclipRealignments = extract(record, false);
        if (softclipRealignments.size() == 0) {
            // nothing to do - just output the record
            coordinateSortedWriter.addAlignment(record);
        } else {
            // perform split read realignment
            SplitReadRealignmentInfo info = new SplitReadRealignmentInfo(record);
            realignments.put(info.alignmentUniqueName, info);
            for (FastqRecord fq : softclipRealignments) {
                aligner.asyncAlign(fq);
                info.outstandingRealignments++;
            }
        }
    }

    private void processAlignmentRecord(SAMRecord supp, Map<String, SplitReadRealignmentInfo> realignments,
                                        SAMFileWriter writer, SAMFileWriter modifiedRecordWriter) throws IOException {
        String lookupKey = SplitReadHelper.getOriginatingAlignmentUniqueName(supp);
        SplitReadRealignmentInfo info = realignments.get(lookupKey);
        if (supp.getSupplementaryAlignmentFlag() || supp.isSecondaryAlignment()) {
            // only consider the best mapping location reported by the aligner
            //log.debug(String.format("%s: ignoring supp alignment", supp.getReadName()));
        } else {
            assert(info.outstandingRealignments > 0);
            info.outstandingRealignments--;
            if (!supp.getReadUnmappedFlag()) {
                info.realignments.add(supp);
                List<FastqRecord> nestedRealignments = extract(supp, true);
                for (FastqRecord fq : nestedRealignments) {
                    aligner.asyncAlign(fq);
                    info.outstandingRealignments++;
                    //log.trace(String.format("%s: performing nested realignment. %d realignments now outstanding", info.originatingRecord.getReadName(), info.outstandingRealignments));
                }
            }
            // all splits identified
            if (info.outstandingRealignments == 0) {
                writeCompletedAlignment(info.originatingRecord, info.realignments, writer, modifiedRecordWriter);
                realignments.remove(lookupKey);
            } else {
                //log.trace(String.format("%s: %d outstanding alignments", info.originatingRecord.getReadName(), info.outstandingRealignments));
            }
        }
    }

    public void createSupplementaryAlignments(final File input, final File output, final File outputModified) throws IOException {
        SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
        SamReaderFactory readerFactory = SamReaderFactory.makeDefault().referenceSequence(pc.getReferenceFile());
        try (SamReader reader = readerFactory.open(input)) {
            boolean unsortedRecordsOutputToSameFile = outputModified == null || output.equals(outputModified);
            SAMFileHeader outputHeader = reader.getFileHeader().clone();
            if (unsortedRecordsOutputToSameFile) {
                outputHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
            }
            SAMFileWriter modifiedWriter = null;
            try (SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(outputHeader, true, output)) {
                if (unsortedRecordsOutputToSameFile) {
                    modifiedWriter = writer;
                } else {
                    SAMFileHeader unsortedHeader = reader.getFileHeader().clone();
                    unsortedHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
                    modifiedWriter = writerFactory.makeSAMOrBAMWriter(SAMFileHeaderUtil.minimal(unsortedHeader), true, outputModified);
                }
                try (AsyncBufferedIterator<SAMRecord> bufferedIt = new AsyncBufferedIterator<>(reader.iterator(), input.getName())) {
                    process(bufferedIt, writer, modifiedWriter);
                }
            } finally {
                CloserUtil.close(modifiedWriter);
            }
        }
    }

    private class SplitReadRealignmentInfo {
        @Override
        public int hashCode() {
            return alignmentUniqueName.hashCode();
        }
        @Override
        public boolean equals(Object obj) {
            return alignmentUniqueName.equals(((SplitReadRealignmentInfo) obj).alignmentUniqueName);
        }
        public SplitReadRealignmentInfo(SAMRecord record) {
            this.originatingRecord = record;
            this.alignmentUniqueName = getEvidenceIdentifierGenerator().getAlignmentUniqueName(record);
        }
        private SAMRecord originatingRecord;
        private String alignmentUniqueName;
        private List<SAMRecord> realignments = new ArrayList<>(2);
        private int outstandingRealignments = 0;
    }
}
