package au.edu.wehi.idsv.alignment;

import au.edu.wehi.idsv.TestHelper;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public abstract class SmithWatermanAlignerTest extends TestHelper {
    protected abstract Aligner create(int match, int mismatch, int ambiguous, int gapOpen, int gapExtend);
    public Aligner create() {
        return create(1, -4, -4, 6, 1);
    }
    @Test
    public void should_align_deletion() {
        Alignment a = create().align_smith_waterman(
                B("AAACCCCCCCCCCCCCTTAATTAATTAATTTT"),
                B("AAACCCCCCCCCCCCCGTTAATTAATTAATTTT"));
        assertEquals("16M1D16M", a.getCigar());
        assertEquals(0, a.getStartPosition());
    }
    @Test
    public void should_use_zero_based_reference_offset() {
        Alignment a = create().align_smith_waterman(B("AACCCTTTTTT"), B("AAACCCTTTTTT"));
        assertEquals("11M", a.getCigar());
        assertEquals(1, a.getStartPosition());
    }
    @Test
    public void should_use_soft_clips() {
        Alignment a = create().align_smith_waterman(
                B("GGGGGGTTTTTT"),
                B( "AAACCCTTTTTT"));
        assertEquals("6S6M", a.getCigar());
        assertEquals(6, a.getStartPosition());
    }
}