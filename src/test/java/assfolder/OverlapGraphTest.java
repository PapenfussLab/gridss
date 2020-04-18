package assfolder;

import au.edu.wehi.idsv.TestHelper;
import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import java.util.List;
import java.util.stream.Collectors;

import static org.junit.Assert.*;

public class OverlapGraphTest extends TestHelper {
    @Test
    public void shouldFindReadOverlaps() {
        List<SAMRecord> records = overlapping(1, 3, 100, 30);
        List<Read> reads = records.stream().map(r -> new Read(r)).collect(Collectors.toList());
        OverlapGraph og = new OverlapGraph(10, 20, 0);
        og.add(reads.get(1));
        og.add(reads.get(2));
        og.add(reads.get(0));
        List<ReadOverlapSuccessor> ros = reads.stream().flatMap(r -> r.overlapSuccessors.stream()).collect(Collectors.toList());
        assertEquals(2, reads.get(0).overlapSuccessors.size());
        assertEquals(1, reads.get(1).overlapSuccessors.size());
        assertEquals(0, reads.get(2).overlapSuccessors.size());
        assertTrue(ros.stream().allMatch(x -> x.mismatches == 0));
    }
    @Test
    public void should_deduplicate_reads() {
        List<SAMRecord> records = overlapping(1, 3, 100, 30);
        Read r = new Read(records.get(0));
        records.get(0).setAlignmentStart(10);
        Read r2 = new Read(records.get(0));
        Read r3 = new Read(records.get(1));
        OverlapGraph og = new OverlapGraph(10, 20, 0);
        assertEquals(r, og.add(r));
        assertEquals(r, og.add(r2));
        assertEquals(r3, og.add(r3));
    }
}