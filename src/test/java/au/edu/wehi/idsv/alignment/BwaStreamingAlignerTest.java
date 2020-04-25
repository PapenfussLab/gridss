package au.edu.wehi.idsv.alignment;

import au.edu.wehi.idsv.TestHelper;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.fastq.FastqRecord;
import org.junit.Test;

import java.io.IOException;
import java.util.Collection;

import static org.junit.Assert.*;

public class BwaStreamingAlignerTest extends TestHelper {
    @Test
    public void should_align_reads() {
        BwaStreamingAligner bwamem = new BwaStreamingAligner(SMALL_FA_FILE, SMALL_FA.getSequenceDictionary(), 2, 10);
        for (FastqRecord fq : ImmutableList.of(
                new FastqRecord("noHit", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
                new FastqRecord("polyA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
                new FastqRecord("polyACGT", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", "", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
                new FastqRecord("random", S(RANDOM).substring(0, 100), "", S(getPolyA(100))))) {
            bwamem.asyncAlign(fq);
        }
        int i = 0;
        while (bwamem.hasAlignmentRecord()) {
            bwamem.getAlignment();
            i++;
        }
        assertTrue(i >= 4);
    }
}