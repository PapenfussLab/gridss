package assfolder;

import au.edu.wehi.idsv.TestHelper;
import com.google.common.collect.Lists;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Collection;

import static org.junit.Assert.*;

public class SeedLookupTest extends TestHelper {
    @Test
    public void should_report_all_matching_positions() {
        Read r1 = new Read(withSequence("AAAAAAAA", Read(0, 1, "8M"))[0]);
        Read r2 = new Read(withSequence("AAAAAAAA", Read(0, 1, "8M"))[0]);
        SeedLookup sl = new SeedLookup(4);
        sl.add(r1);
        Collection<ReadOffset> col = sl.findOverlaps(r2);
        //     AAAAAAAA
        // AAAAAAAA
        //  AAAAAAAA
        //   AAAAAAAA
        //    AAAAAAAA
        //     AAAAAAAA
        //      AAAAAAAA
        //       AAAAAAAA
        //        AAAAAAAA
        //         AAAAAAAA
        assertEquals(9, col.size());
        for (int i = -4; i <= 4 ; i++) {
            assertTrue(col.contains(new ReadOffset(r1, i)));
        }
    }
    @Test
    public void should_report_offset_relative_to_read() {
        Read r1 = new Read(withSequence("ATTAAATC", Read(0, 1, "8M"))[0]);
        Read r2 = new Read(withSequence("CCAT", Read(0, 1, "4M"))[0]);
        SeedLookup sl = new SeedLookup(2);
        sl.add(r1);
        Collection<ReadOffset> col = sl.findOverlaps(r2);
        assertEquals(2, col.size());
        //   01234567
        //   ATTAAATC
        // CCAT
        //      CCAT
        assertTrue(col.contains(new ReadOffset(r1, 2)));
        assertTrue(col.contains(new ReadOffset(r1, -3)));
    }
    @Test
    public void offset_should_be_number_of_bases_result_read_is_after_arg_read() {
        Read r1 = new Read(withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(2, 1, "100M"))[0]);
        Read r2 = new Read(withSequence( "ATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGATT", Read(2, 2, "100M"))[0]);
        SeedLookup sl = new SeedLookup(16);
        sl.add(r1);
        sl.add(r2);
        ArrayList<ReadOffset> list = Lists.newArrayList(sl.findOverlaps(r1));
        assertEquals(1, list.size());
        assertEquals(1, list.get(0).offset);
        assertEquals(r2, list.get(0).read);
        list = Lists.newArrayList(sl.findOverlaps(r2));
        assertEquals(1, list.size());
        assertEquals(-1, list.get(0).offset);
        assertEquals(r1, list.get(0).read);
    }
}
