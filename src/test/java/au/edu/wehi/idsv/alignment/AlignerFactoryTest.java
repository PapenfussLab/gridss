package au.edu.wehi.idsv.alignment;

import au.edu.wehi.idsv.LinuxTests;
import org.junit.Test;
import org.junit.experimental.categories.Category;

import static org.junit.Assert.assertFalse;

public class AlignerFactoryTest {
    /**
     * Intel's aligner is broken and unusable:
     * https://github.com/Intel-HLS/GKL/issues/104
     */
    @Test
    @Category(LinuxTests.class)
    public void should_not_default_to_GKL_aligner() {
        Aligner aligner = AlignerFactory.create();
        assertFalse(aligner instanceof GKLAligner);
    }
}