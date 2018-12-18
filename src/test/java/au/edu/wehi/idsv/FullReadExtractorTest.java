package au.edu.wehi.idsv;

import au.edu.wehi.idsv.bed.IntervalBed;
import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.*;

public class FullReadExtractorTest extends TestHelper {
    private static LinearGenomicCoordinate lgc = new PaddedLinearGenomicCoordinate(getSequenceDictionary(), LCCB);
    public class FullReadExtractorStub extends FullReadExtractor {
        public FullReadExtractorStub(LinearGenomicCoordinate lgc, IntervalBed bed, boolean mates, boolean splits) {
            super(lgc, bed, mates, splits);
        }

        @Override
        public void extract(File input, File output) throws IOException {

        }
    }
    @Test
    public void shouldExtractOverlappingRecords() {
        IntervalBed bed = new IntervalBed(lgc.getDictionary(), lgc);
        bed.addInterval(0, 10, 20);
        bed.addInterval(0, 30, 40);
        FullReadExtractorStub fre = new FullReadExtractorStub(lgc, bed, false, false);
        assertFalse(fre.shouldExtract(Read(0, 1, "1M")));
        assertFalse(fre.shouldExtract(Read(0, 9, "1M")));
        assertTrue(fre.shouldExtract(Read(0, 9, "2M")));
        assertTrue(fre.shouldExtract(Read(0, 10, "1M")));
        assertTrue(fre.shouldExtract(Read(0, 20, "1M")));
        assertFalse(fre.shouldExtract(Read(0, 21, "1M")));
        assertTrue(fre.shouldExtract(Read(0, 35, "1M")));
        assertFalse(fre.shouldExtract(Read(0, 45, "1M")));
    }
    @Test
    public void should_extract_if_mate_overlapping() {
        IntervalBed bed = new IntervalBed(lgc.getDictionary(), lgc);
        bed.addInterval(0, 10, 20);
        bed.addInterval(0, 30, 40);
        FullReadExtractorStub fre = new FullReadExtractorStub(lgc, bed, true, false);
        assertFalse(fre.shouldExtract(DP(1, 1, "100M", true, 0, 1, "1M", true)[0]));
        assertFalse(fre.shouldExtract(DP(1, 1, "100M", true, 0, 9, "1M", true)[0]));
        assertTrue(fre.shouldExtract(DP(1, 1, "100M", true, 0, 9, "2M", true)[0]));
        assertTrue(fre.shouldExtract(DP(1, 1, "100M", true, 0, 10, "1M", true)[0]));
        assertTrue(fre.shouldExtract(DP(1, 1, "100M", true, 0, 20, "1M", true)[0]));
        assertFalse(fre.shouldExtract(DP(1, 1, "100M", true, 0, 21, "1M", true)[0]));
        assertTrue(fre.shouldExtract(DP(1, 1, "100M", true, 0, 35, "1M", true)[0]));
        assertFalse(fre.shouldExtract(DP(1, 1, "100M", true, 0, 45, "1M", true)[0]));

        fre = new FullReadExtractorStub(lgc, bed, false, true);
        assertFalse(fre.shouldExtract(DP(1, 1, "100M", true, 0, 9, "2M", true)[0]));
    }
    @Test
    public void should_extract_if_split_overlapping() {
        IntervalBed bed = new IntervalBed(lgc.getDictionary(), lgc);
        bed.addInterval(0, 100, 200);
        bed.addInterval(1, 30, 40);
        FullReadExtractorStub fre = new FullReadExtractorStub(lgc, bed, false, true);
        //(rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+
        assertTrue(fre.shouldExtract(withAttr("SA", "polyA,95,+,5S10M,0,0",Read(0, 1, "5M10S"))[0]));
        assertFalse(fre.shouldExtract(withAttr("SA", "polyA,50,+,5S10M,0,0",Read(0, 1, "5M10S"))[0]));
        assertTrue(fre.shouldExtract(withAttr("SA", "polyA,50,+,5S10M,0,0;polyA,95,+,5S10M,0,0",Read(0, 1, "5M10S"))[0]));

        fre = new FullReadExtractorStub(lgc, bed, true, false);
        assertFalse(fre.shouldExtract(withAttr("SA", "polyA,50,+,5S10M,0,0;polyA,95,+,5S10M,0,0",Read(0, 1, "5M10S"))[0]));
    }
}