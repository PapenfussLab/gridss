package au.edu.wehi.idsv.alignment;

import com.intel.gkl.smithwaterman.IntelSmithWaterman;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWNativeAlignerResult;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.junit.Ignore;
import org.junit.Test;

import static org.junit.Assert.*;

@Ignore("It's broken: https://github.com/Intel-HLS/GKL/issues/104")
public class GKLAlignerTest extends SmithWatermanAlignerTest {
    @Override
    protected Aligner create(int match, int mismatch, int ambiguous, int gapOpen, int gapExtend) {
        IntelSmithWaterman isw = new IntelSmithWaterman();
        if (!isw.load(null)) {
            throw new RuntimeException("Unable to load GKL native library");
        }
        return new GKLAligner(match, mismatch, gapOpen, gapExtend, isw);
    }
    @Test
    public void should_return_maximal_alignment_match() {
        IntelSmithWaterman isw = new IntelSmithWaterman();
        assertTrue(isw.load(null));
        SWNativeAlignerResult alignment = isw.align(
                "GGGGGGTTTTTT".getBytes(),
                "AAACCCTTTTTT".getBytes(),
                // BWA alignment scoring matrix
                new SWParameters(1, -4, -6, -1),
                SWOverhangStrategy.SOFTCLIP);
        assertEquals("6S6M", alignment.cigar);
        assertEquals(6, alignment.alignment_offset);
    }
}