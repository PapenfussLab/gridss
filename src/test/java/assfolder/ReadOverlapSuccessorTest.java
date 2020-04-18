package assfolder;

import au.edu.wehi.idsv.TestHelper;
import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import static org.junit.Assert.*;

public class ReadOverlapSuccessorTest extends TestHelper {
    @Test
    public void should_report_mismatches() {
        List<SAMRecord> records = overlapping(1, 2, 100, 30);
        List<Read> reads = records.stream().map(r -> new Read(r)).collect(Collectors.toList());
        ReadOverlapSuccessor ros = new ReadOverlapSuccessor(reads.get(0), reads.get(1), 30);
        assertEquals(70, ros.overlapLength);
        assertEquals(0, ros.mismatches);

        ros = new ReadOverlapSuccessor(reads.get(1), reads.get(0), 30);
        assertEquals(70, ros.overlapLength);
    }
    @Test
    public void should_report_contained_reads_with_minus_sequence_prefix() {
        Read r1 = new Read(withSequence(   "ATTAAATC", Read(0, 1, "8M"))[0]);
        Read r2 = new Read(withSequence("CCCATTAAATCTT", Read(0, 1, "13M"))[0]);
        ReadOverlapSuccessor ros = new ReadOverlapSuccessor(r2, r1, 3);
        assertEquals(3, ros.getLeadingBases());
        assertEquals(-2, ros.getTrailingBases());
        assertEquals("CCC", ros.getPrefixSequence(r2));
        assertEquals("-TT",ros.getSuffixSequence(r2));
    }
    @Test
    public void should_report_dangling_sequence() {
        Read r1 = new Read(withSequence(   "ATTAAATC", Read(0, 1, "8M"))[0]);
        Read r2 = new Read(withSequence("CCCATT", Read(0, 1, "6M"))[0]);
        ReadOverlapSuccessor ros = new ReadOverlapSuccessor(r2, r1, 3);
        assertEquals(3, ros.getLeadingBases());
        assertEquals(5, ros.getTrailingBases());
        assertEquals("CCC", ros.getPrefixSequence(r2));
        assertEquals("AAATC",ros.getSuffixSequence(r2));
    }
}