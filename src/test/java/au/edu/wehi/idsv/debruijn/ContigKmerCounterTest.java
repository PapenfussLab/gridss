package au.edu.wehi.idsv.debruijn;

import au.edu.wehi.idsv.TestHelper;
import com.google.common.collect.ImmutableList;
import org.junit.Assert;
import org.junit.Test;

public class ContigKmerCounterTest extends TestHelper {
    @Test
    public void should_count_at_end_of_read_and_reference() {
        ContigKmerCounter ckc = new ContigKmerCounter(ImmutableList.of("test1", "test2"), ImmutableList.of(B("AACCGGTTT"), B("A")), 4, 1);
        Assert.assertEquals(1, ckc.count(B("CGTTT")));
        Assert.assertEquals("test1", ckc.getContigs().get(0));
        Assert.assertEquals("test2", ckc.getContigs().get(1));
        Assert.assertEquals(1, ckc.getKmerCounts().getLong(0));
        Assert.assertEquals(0, ckc.getKmerCounts().getLong(1));
    }
    @Test
    public void should_ignore_duplicate_contigs() {
        ContigKmerCounter ckc = new ContigKmerCounter(ImmutableList.of("test1", "test2", "test1"), ImmutableList.of(B("AACCGGTTT"), B("A"), B("AAAAAAAAAAAAAAAAAA")), 3, 1);
        Assert.assertEquals(2, ckc.getContigs().size());
        Assert.assertEquals(1, ckc.count(B("AAA")));
        Assert.assertEquals("test1", ckc.getContigs().get(0));
        Assert.assertEquals("test2", ckc.getContigs().get(1));
    }
    @Test
    public void should_not_double_count_repeated_kmers() {
        ContigKmerCounter ckc = new ContigKmerCounter(ImmutableList.of("test1", "test2"), ImmutableList.of(B("CCTCCTT"), B("A")), 3, 1);
        Assert.assertEquals(3, ckc.count(B("CTCCT")));
    }
    @Test
    public void should_incorporate_stride() {
        ContigKmerCounter ckc = new ContigKmerCounter(ImmutableList.of("test1", "test2"),
            ImmutableList.of(
                B("CCACCCA"),
                B("AACCA")),
            3, 3);
        Assert.assertEquals(1, ckc.count(B("CCA")));
        Assert.assertEquals(1, ckc.count(B("AAC")));
        Assert.assertEquals(0, ckc.count(B("CAC")));
        Assert.assertEquals(1, ckc.count(B("CCC")));
        Assert.assertEquals(2, ckc.getKmerCounts().getLong(0));
        Assert.assertEquals(1, ckc.getKmerCounts().getLong(1));
    }
}