package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.SingleReadEvidence;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.positional.KmerEvidence;
import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import java.util.List;
import java.util.stream.Collectors;

import static org.junit.Assert.*;

public class KmerEvidenceSetTest extends TestHelper {
    @Test
    public void test_operations() {
        MockSAMEvidenceSource ses = SES();
        SAMRecord r1 = Read(0, 1, "10M5D10M15S");
        int k = 2;
        List<KmerEvidence> list = SingleReadEvidence.createEvidence(SES(), 1, r1).stream()
                .map(sre -> KmerEvidence.create(k, sre))
                .collect(Collectors.toList());
        KmerEvidenceSet s = new KmerEvidenceSet();
        assertEquals(0, s.size());
        assertTrue(s.add(list.get(0)));
        assertEquals(1, s.size());
        assertTrue(s.contains(list.get(0)));
        assertFalse(s.contains(list.get(1)));
        assertFalse(s.add(list.get(0)));
        assertTrue(s.remove(list.get(0)));
        assertEquals(0, s.size());
    }
}