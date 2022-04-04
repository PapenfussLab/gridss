package au.edu.wehi.idsv.sim;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.vcf.SvType;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


public class SimpleEventTest extends TestHelper {
    @Test
    public void should_encode_dup() {
        SimpleEvent e = new SimpleEvent(SvType.DUP, 1, 3, 1, "");
        assertEquals(3, e.getGenomicWidth());
        assertEquals("CGTCGT", e.getVariantSeq(SMALL_FA, 0, 0));
        assertEquals("ACGTCGTAC", e.getVariantSeq(SMALL_FA, 1, 2));
        assertEquals("A", e.asVariantContextBuilder(SMALL_FA, true).make().getReference().getBaseString());
        assertEquals("<DUP>", e.asVariantContextBuilder(SMALL_FA, true).make().getAlternateAllele(0).getDisplayString());
        assertEquals("ACGT", e.asVariantContextBuilder(SMALL_FA, false).make().getReference().getBaseString());
        assertEquals("ACGTCGT", e.asVariantContextBuilder(SMALL_FA, false).make().getAlternateAllele(0).getBaseString());

    }
    @Test
    public void should_encode_ins() {
        SimpleEvent e = new SimpleEvent(SvType.INS, 1, 3, 1, "TTT");
        assertEquals(0, e.getGenomicWidth());
        assertEquals("TTT", e.getVariantSeq(SMALL_FA, 0, 0));
        assertEquals("ATTTCG", e.getVariantSeq(SMALL_FA, 1, 2));
        assertEquals("A", e.asVariantContextBuilder(SMALL_FA, true).make().getReference().getBaseString());
        assertEquals("<INS>", e.asVariantContextBuilder(SMALL_FA, true).make().getAlternateAllele(0).getDisplayString());
        assertEquals("A", e.asVariantContextBuilder(SMALL_FA, false).make().getReference().getBaseString());
        assertEquals("ATTT", e.asVariantContextBuilder(SMALL_FA, false).make().getAlternateAllele(0).getBaseString());
    }
    @Test
    public void should_encode_inv() {
        SimpleEvent e = new SimpleEvent(SvType.INV, 1, 3, 1, "");
        assertEquals(3, e.getGenomicWidth());
        assertEquals("ACG", e.getVariantSeq(SMALL_FA, 0, 0));
        assertEquals("AACGAC", e.getVariantSeq(SMALL_FA, 1, 2));
        assertEquals("A", e.asVariantContextBuilder(SMALL_FA, true).make().getReference().getBaseString());
        assertEquals("<INV>", e.asVariantContextBuilder(SMALL_FA, true).make().getAlternateAllele(0).getDisplayString());
        assertEquals("ACGT", e.asVariantContextBuilder(SMALL_FA, false).make().getReference().getBaseString());
        assertEquals("AACG", e.asVariantContextBuilder(SMALL_FA, false).make().getAlternateAllele(0).getBaseString());
    }
    @Test
    public void should_encode_del() {
        SimpleEvent e = new SimpleEvent(SvType.DEL, 1, 3, 1, "");
        assertEquals(3, e.getGenomicWidth());
        assertEquals("", e.getVariantSeq(SMALL_FA, 0, 0));
        assertEquals("AAC", e.getVariantSeq(SMALL_FA, 1, 2));
        assertEquals("A", e.asVariantContextBuilder(SMALL_FA, true).make().getReference().getBaseString());
        assertEquals("<DEL>", e.asVariantContextBuilder(SMALL_FA, true).make().getAlternateAllele(0).getDisplayString());
        assertEquals("ACGT", e.asVariantContextBuilder(SMALL_FA, false).make().getReference().getBaseString());
        assertEquals("A", e.asVariantContextBuilder(SMALL_FA, false).make().getAlternateAllele(0).getBaseString());
    }
    @Test
    public void should_flag_reference_sequence() {
        SimpleEvent e = new SimpleEvent(SvType.INV, 1, 2, 1, "");
        assertTrue(e.asVariantContextBuilder(SMALL_FA, false).getFilters().contains("REF"));
    }
}