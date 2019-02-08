package au.edu.wehi.idsv;

import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.sam.SAMFileUtil;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import static au.edu.wehi.idsv.sam.ChimericAlignment.getChimericAlignments;

public class IndexedReadExtractor extends ReadExtractor {
    private static final Log log = Log.getInstance(IndexedReadExtractor.class);
    public IndexedReadExtractor(LinearGenomicCoordinate lgc, IntervalBed bed, boolean extractMates, boolean extractSplits) {
        super(lgc, bed, extractMates, extractSplits);
    }

    @Override
    public void extract(File input, File output, int workerThreads) throws IOException {
        extract(input, output, getRegionBed().asQueryInterval());
    }
    private void extract(File input, File output, QueryInterval[] intervals) throws IOException {
        IntervalBed remoteLocations = new IntervalBed(getLinearGenomicCoordinate());
        boolean shouldLookupUnmapped = false;
        File regionOut = FileSystemContext.getWorkingFileFor(output);
        File offTargetOut =  FileSystemContext.getWorkingFileFor(output, "mate_splits");
        try (SamReader reader = SamReaderFactory.makeDefault().open(input)) {
            if (!reader.hasIndex()) {
                throw new RuntimeException("Missing BAM index for " + input.getName());
            }
            SAMFileHeader header = reader.getFileHeader();
            log.info(String.format("Extracting %d intervals.", intervals.length));
            try (SAMRecordIterator it = reader.query(intervals, false)) {
                //ProgressLoggingSAMRecordIterator it = new ProgressLoggingSAMRecordIterator(it, new ProgressLogger(log));
                try (SAMFileWriter writer = new SAMFileWriterFactory().setCompressionLevel(0).makeBAMWriter(header, true, regionOut)) {
                    while (it.hasNext()) {
                        SAMRecord r = it.next();
                        if (overlapsRegionBed(r)) {
                            writer.addAlignment(r);
                            if (shouldExtractMates() && r.getReadPairedFlag()) {
                                if (r.getMateUnmappedFlag()) {
                                    shouldLookupUnmapped = true;
                                } else {
                                    remoteLocations.addInterval(r.getMateReferenceIndex(), r.getMateAlignmentStart(), r.getMateAlignmentStart());
                                }
                            }
                            List<ChimericAlignment> splits = getChimericAlignments(r);
                            if (shouldExtractSplits() && !splits.isEmpty()) {
                                for (ChimericAlignment ca : splits) {
                                    remoteLocations.addInterval(getLinearGenomicCoordinate().getDictionary().getSequenceIndex(ca.rname), ca.pos, ca.pos);
                                }
                            }
                        }
                    }
                }
            }
            // iterator over remote targets
            QueryInterval[] offTarget = remoteLocations.asQueryInterval();
            log.info(String.format("Querying %d intervals for mates and split reads.", offTarget.length));
            try (SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, offTargetOut)) {
                try (SAMRecordIterator it = reader.query(remoteLocations.asQueryInterval(), false)) {
                    while (it.hasNext()) {
                        SAMRecord r = it.next();
                        if (!overlapsRegionBed(r) && shouldExtract(r)) {
                            writer.addAlignment(r);
                        }
                    }
                }
                if (shouldLookupUnmapped) {
                    try (SAMRecordIterator it = reader.queryUnmapped()) {
                        while (it.hasNext()) {
                            SAMRecord r = it.next();
                            if (shouldExtract(r)) {
                                writer.addAlignment(r);
                            }
                        }
                    }
                }
            }
        }
        SAMFileUtil.merge(ImmutableList.of(regionOut, offTargetOut), output);
        Files.delete(regionOut.toPath());
        Files.delete(offTargetOut.toPath());
    }
}
