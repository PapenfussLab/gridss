package au.edu.wehi.idsv.sim;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


import java.nio.charset.StandardCharsets;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class RandomReferenceSequenceGeneratorTest {
    @Test
    public void should_sample_from_reference() {
        byte[] ref = "ACGTGTGA".getBytes(StandardCharsets.US_ASCII);
        RandomReferenceSequenceGenerator gen = new RandomReferenceSequenceGenerator(0, ref);
        Set<String> results = IntStream.range(0, 128).mapToObj(x -> new String(gen.getBases(6))).collect(Collectors.toSet());
        assertEquals(3, results.size());
        assertTrue(results.contains("ACGTGT"));
        assertTrue(results.contains("CGTGTG"));
        assertTrue(results.contains("GTGTGA"));
    }
}