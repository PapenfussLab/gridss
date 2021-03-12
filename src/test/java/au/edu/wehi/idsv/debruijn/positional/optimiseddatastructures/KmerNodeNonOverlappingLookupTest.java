package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.positional.ImmutableKmerNode;
import au.edu.wehi.idsv.debruijn.positional.KmerNode;
import au.edu.wehi.idsv.debruijn.positional.KmerPathNode;
import com.google.common.collect.ImmutableList;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

public class KmerNodeNonOverlappingLookupTest extends TestHelper {
    @Test
    public void getUnique_should_return_matching() {
        int k = 4;
        List<ImmutableKmerNode> nodes = ImmutableList.of(
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 1, 1, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AATT")), 2, 2, false, 1));
        KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
        nodes.stream().forEach(n -> lookup.add(n));
        lookup.sanityCheck();
        Assert.assertEquals(nodes.get(0), lookup.getUniquePredecessor(nodes.get(1)));
        Assert.assertNull(lookup.getUniquePredecessor(nodes.get(0)));
        Assert.assertEquals(nodes.get(1), lookup.getUniqueSuccessor(nodes.get(0)));
        Assert.assertNull(lookup.getUniqueSuccessor(nodes.get(1)));
    }
    @Test
    public void getUniquePredecessor_should_return_matching() {
        int k = 4;
        List<ImmutableKmerNode> nodes = ImmutableList.of(
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 1, 1, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AATT")), 2, 2, false, 1));
        KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
        nodes.stream().forEach(n -> lookup.add(n));
        lookup.sanityCheck();
        Assert.assertEquals(nodes.get(0), lookup.getUniquePredecessor(nodes.get(1)));
        Assert.assertNull(lookup.getUniquePredecessor(nodes.get(0)));
    }
    @Test
    public void getUniquePredecessor_should_not_return_mismatched() {
        int k = 4;
        // contained
        List<ImmutableKmerNode> nodes = ImmutableList.of(
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 1, 1, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AATT")), 0, 3, false, 1));
        {
            final KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
            nodes.stream().forEach(n -> lookup.add(n));
            lookup.sanityCheck();
            Assert.assertFalse(nodes.stream().anyMatch(n -> lookup.getUniquePredecessor(n) != null));
        }
        // containing
        nodes = ImmutableList.of(
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 1, 4, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AATT")), 2, 3, false, 1));
        {
            final KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
            nodes.stream().forEach(n -> lookup.add(n));
            lookup.sanityCheck();
            Assert.assertFalse(nodes.stream().anyMatch(n -> lookup.getUniquePredecessor(n) != null));
        }
        // before
        nodes = ImmutableList.of(
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 1, 1, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AATT")), 1, 1, false, 1));
        {
            final KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
            nodes.stream().forEach(n -> lookup.add(n));
            lookup.sanityCheck();
            Assert.assertFalse(nodes.stream().anyMatch(n -> lookup.getUniquePredecessor(n) != null));
        }
        // after
        nodes = ImmutableList.of(
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 1, 1, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AATT")), 1, 1, false, 1));
        {
            final KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
            nodes.stream().forEach(n -> lookup.add(n));
            lookup.sanityCheck();
            Assert.assertFalse(nodes.stream().anyMatch(n -> lookup.getUniquePredecessor(n) != null));
        }
        // multiple position
        nodes = ImmutableList.of(
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 1, 1, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 2, 2, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AATT")), 2, 3, false, 1));
        {
            final KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
            nodes.stream().forEach(n -> lookup.add(n));
            lookup.sanityCheck();
            Assert.assertFalse(nodes.stream().anyMatch(n -> lookup.getUniquePredecessor(n) != null));
        }
        // multiple kmer
        nodes = ImmutableList.of(
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 1, 1, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("CAAT")), 1, 1, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AATT")), 2, 2, false, 1));
        {
            final KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
            nodes.stream().forEach(n -> lookup.add(n));
            lookup.sanityCheck();
            Assert.assertFalse(nodes.stream().anyMatch(n -> lookup.getUniquePredecessor(n) != null));
        }
    }
    @Test
    public void next_should_match_prev() {
        List<ImmutableKmerNode> nodes = new ArrayList<>();
        int k = 3;
        for (long kmer  = 0; kmer < 1 << (2 * k); kmer++) {
            nodes.add(new ImmutableKmerNode(kmer, 1, 2, false, 1));
            nodes.add(new ImmutableKmerNode(kmer, 3, 3, false, 1));
            nodes.add(new ImmutableKmerNode(kmer, 4, 4, false, 1));
            nodes.add(new ImmutableKmerNode(kmer, 5, 10, false, 1));
        }
        KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
        nodes.stream().forEach(n -> lookup.add(n));
        lookup.sanityCheck();
        for (ImmutableKmerNode n : nodes) {
            for (KmerNode next : lookup.nextNodes(n)) {
                Assert.assertTrue(lookup.prevNodes(next).contains(n));
            }
        }
    }
    @Test
    public void next_should_include_matching_nodes() {
        int k = 4;
        List<ImmutableKmerNode> nodes = ImmutableList.of(
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAA")), 1, 3, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 1, 1, false, 1), // no overlap
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 2, 2, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 4, 10, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAC")), 1, 2, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAG")), 2, 4, false, 1));
        final KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
        nodes.stream().forEach(n -> lookup.add(n));
        lookup.sanityCheck();
        List<KmerNode> x = lookup.nextNodes(nodes.get(0));
        Assert.assertEquals(5, x.size());
    }
    @Test
    public void adjustForMerge_should_remove_internal_lookups() {
        int k = 4;
        List<ImmutableKmerNode> nodes = ImmutableList.of(
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("TAAA")), 1, 1, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 2, 2, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AATT")), 3, 3, false, 1));
        final KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
        nodes.stream().forEach(n -> lookup.add(n));
        lookup.sanityCheck();
        lookup.remove(nodes.get(0));
        lookup.sanityCheck();
        KmerPathNode kpn = new KmerPathNode(nodes.get(0).firstKmer(), 1, 1, false, 1);
        lookup.add(kpn);
        lookup.sanityCheck();
        for (int i = 1; i <= 2 ; i++) {
            lookup.adjustForMerge(kpn, nodes.get(i));
            // sanity check will fail here
            kpn.append(nodes.get(i));
            // and now should be good again
            lookup.sanityCheck();
        }
    }
    @Test
    public void next_should_include_end_overlap() {
        int k = 4;
        List<ImmutableKmerNode> nodes = ImmutableList.of(
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("GAAA")), 1, 3, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 4, 10, false, 1));
        final KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
        nodes.stream().forEach(n -> lookup.add(n));
        lookup.sanityCheck();
        List<KmerNode> x = lookup.nextNodes(nodes.get(0));
        Assert.assertEquals(1, x.size());
    }
}