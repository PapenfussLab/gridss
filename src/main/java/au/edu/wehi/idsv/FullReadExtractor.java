package au.edu.wehi.idsv;

import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.FileHelper;
import com.google.common.collect.Range;
import gridss.ExtractFullReads;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static au.edu.wehi.idsv.sam.ChimericAlignment.getChimericAlignments;
import static htsjdk.samtools.SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX;
import static htsjdk.samtools.SAMRecord.NO_ALIGNMENT_START;

public class FullReadExtractor {
    private static final Log log = Log.getInstance(FullReadExtractor.class);
    private final LinearGenomicCoordinate lgc;
    private final IntervalBed bed;
    private final boolean extractMates;
    private final boolean extractSplits;

    public FullReadExtractor(LinearGenomicCoordinate lgc, IntervalBed bed, boolean extractMates, boolean extractSplits) {
        this.lgc = lgc;
        this.bed = bed;
        this.extractMates = extractMates;
        this.extractSplits = extractSplits;
    }

    public LinearGenomicCoordinate getLinearGenomicCoordinate() {
        return lgc;
    }

    public IntervalBed getRegionBed() {
        return bed;
    }

    public boolean overlapsRegionBed(SAMRecord r) {
        return overlapsRegionBed(r, 0);
    }

    public boolean overlapsRegionBed(SAMRecord r, int padding) {
        return !r.getReadUnmappedFlag() && bed.overlaps(r.getReferenceIndex(), r.getStart() - padding, r.getEnd() + padding);
    }
    protected Range<Long> getRange(SAMRecord r) {
        int refIndex = r.getReferenceIndex();
        int start = r.getAlignmentStart();
        if (refIndex == NO_ALIGNMENT_REFERENCE_INDEX || start == NO_ALIGNMENT_START) {
            return null;
        } else {
            return rangeOf(lgc.getLinearCoordinate(refIndex, start), r.getCigar());
        }
    }
    protected Range<Long> getMateRange(SAMRecord r) {
        if (!r.getReadPairedFlag() || r.getMateUnmappedFlag()) {
            return null;
        }
        Object mc = r.getAttribute(SAMTag.MC.name());
        Cigar mateCigar = null;
        if (mc instanceof String) {
            mateCigar = TextCigarCodec.decode((String)mc);
        }
        long start = lgc.getLinearCoordinate(r.getMateReferenceIndex(), r.getMateAlignmentStart());
        return rangeOf(start, mateCigar);
    }

    protected List<Range<Long>> getSplitRanges(SAMRecord r) {
        List<ChimericAlignment> cal = getChimericAlignments(r);
        List<Range<Long>> result = new ArrayList<>(cal.size());
        for (ChimericAlignment ca : cal) {
            long startOffset = lgc.getLinearCoordinate(ca.rname, ca.pos);
            result.add(rangeOf(startOffset, ca.cigar));
        }
        return result;
    }
    protected Range<Long> rangeOf(long startOffset, Cigar cigar) {
        int refLength = (cigar == null || cigar.isEmpty()) ? 1 : cigar.getReferenceLength();
        return Range.closedOpen(startOffset, startOffset + refLength);
    }

    public boolean shouldExtract(SAMRecord r) {
        if (overlapsRegionBed(r)) {
            return true;
        }
        if (shouldExtractMates()) {
            Range<Long> mate = getMateRange(r);
            if (mate != null && bed.overlaps(mate)) {
                return true;
            }
        }
        if (shouldExtractSplits()) {
            for (Range<Long> chim : getSplitRanges(r)) {
                if (bed.overlaps(chim)) {
                    return true;
                }
            }
        }
        return false;
    }

    public boolean shouldExtractMates() {
        return extractMates;
    }

    public boolean shouldExtractSplits() {
        return extractSplits;
    }

    public void extract(File input, File output, int workerThreads) throws IOException {
        File tmpOut = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(output) : output;
        try (SamReader reader = SamReaderFactory.makeDefault().open(input)) {
            SAMFileHeader header = reader.getFileHeader();
            try (AsyncBufferedIterator<SAMRecord> asyncIt = new AsyncBufferedIterator<>(reader.iterator(), input.getName())) {
                ProgressLoggingSAMRecordIterator it = new ProgressLoggingSAMRecordIterator(asyncIt, new ProgressLogger(log));
                try (SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, tmpOut)) {
                    while (it.hasNext()) {
                        SAMRecord r = it.next();
                        if (shouldExtract(r)) {
                            writer.addAlignment(r);
                        }
                    }
                }
            }
        }
        if (tmpOut != output) {
            FileHelper.move(tmpOut, output, true);
        }
    }
}
