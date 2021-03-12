package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.positional.ImmutableKmerNode;
import au.edu.wehi.idsv.debruijn.positional.KmerNode;
import au.edu.wehi.idsv.debruijn.positional.KmerNodeUtil;
import au.edu.wehi.idsv.debruijn.positional.KmerPathNode;
import au.edu.wehi.idsv.util.IntervalUtil;
import com.google.common.collect.ImmutableList;
import org.junit.Assert;
import org.junit.Test;

import java.util.*;
import java.util.stream.Collectors;

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
        Assert.assertEquals(nodes.get(0), lookup.getUniqueFullWidthPredecessor(nodes.get(1)));
        Assert.assertNull(lookup.getUniqueFullWidthPredecessor(nodes.get(0)));
        Assert.assertEquals(nodes.get(1), lookup.getUniqueFullWidthSuccessor(nodes.get(0)));
        Assert.assertNull(lookup.getUniqueFullWidthSuccessor(nodes.get(1)));
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
        Assert.assertEquals(nodes.get(0), lookup.getUniqueFullWidthPredecessor(nodes.get(1)));
        Assert.assertNull(lookup.getUniqueFullWidthPredecessor(nodes.get(0)));
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
            Assert.assertFalse(nodes.stream().anyMatch(n -> lookup.getUniqueFullWidthPredecessor(n) != null));
        }
        // containing
        nodes = ImmutableList.of(
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 1, 4, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AATT")), 2, 3, false, 1));
        {
            final KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
            nodes.stream().forEach(n -> lookup.add(n));
            lookup.sanityCheck();
            Assert.assertFalse(nodes.stream().anyMatch(n -> lookup.getUniqueFullWidthPredecessor(n) != null));
        }
        // before
        nodes = ImmutableList.of(
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 1, 1, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AATT")), 1, 1, false, 1));
        {
            final KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
            nodes.stream().forEach(n -> lookup.add(n));
            lookup.sanityCheck();
            Assert.assertFalse(nodes.stream().anyMatch(n -> lookup.getUniqueFullWidthPredecessor(n) != null));
        }
        // after
        nodes = ImmutableList.of(
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 1, 1, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AATT")), 1, 1, false, 1));
        {
            final KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
            nodes.stream().forEach(n -> lookup.add(n));
            lookup.sanityCheck();
            Assert.assertFalse(nodes.stream().anyMatch(n -> lookup.getUniqueFullWidthPredecessor(n) != null));
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
            Assert.assertFalse(nodes.stream().anyMatch(n -> lookup.getUniqueFullWidthPredecessor(n) != null));
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
            Assert.assertFalse(nodes.stream().anyMatch(n -> lookup.getUniqueFullWidthPredecessor(n) != null));
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
                List<KmerNode> prevOfNext = lookup.prevNodes(next);
                Assert.assertTrue(prevOfNext.contains(n));
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
    @Test
    public void replace_should_remove_old_node_from_lookup() {
        int k = 4;
        List<ImmutableKmerNode> nodes = ImmutableList.of(
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("GAAA")), 1, 3, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 4, 10, false, 1));
        final KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
        nodes.stream().forEach(n -> lookup.add(n));
        lookup.sanityCheck();
        KmerPathNode kpn = new KmerPathNode(nodes.get(0).firstKmer(), 1, 3, false, 1);
        lookup.replace(nodes.get(0), kpn);
        lookup.sanityCheck();
        Assert.assertTrue(kpn == lookup.prevNodes(nodes.get(1)).get(0));
    }
    @Test
    public void should_return_successor() {
        int k = 4;
        List<ImmutableKmerNode> nodes = ImmutableList.of(
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("TAAA")), 1, 1, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAA")), 2, 2, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAA")), 3, 3, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 4, 4, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAC")), 5, 5, false, 1));
        final KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
        nodes.stream().forEach(n -> lookup.add(n));
        lookup.sanityCheck();
        Assert.assertTrue(nodes.get(0) == lookup.getUniqueFullWidthPredecessor(nodes.get(1)));
        Assert.assertTrue(nodes.get(1) == lookup.getUniqueFullWidthSuccessor(nodes.get(0)));
    }
    @Test
    public void uniques_should_match() {
        int k = 4;
        List<ImmutableKmerNode> nodes = ImmutableList.of(
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("TAAG")), 1, 5, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAGA")), 2, 6, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAG")), 0, 0, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAA")), 7, 7, false, 1));
        final KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
        nodes.stream().forEach(n -> lookup.add(n));
        lookup.sanityCheck();
        Assert.assertTrue(nodes.get(0) == lookup.getUniqueFullWidthPredecessor(nodes.get(1)));
        Assert.assertTrue(nodes.get(1) == lookup.getUniqueFullWidthSuccessor(nodes.get(0)));
    }
    @Test
    public void nextNodes_bounds_should_match() {
        int k = 4;
        List<ImmutableKmerNode> nodes = ImmutableList.of(
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("TAAA")), 1, 4, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAA")), 0, 0, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAA")), 1, 1, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAA")), 2, 2, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAA")), 3, 3, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAA")), 4, 4, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAA")), 5, 5, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAA")), 6, 6, false, 1));
        final KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
        nodes.stream().forEach(n -> lookup.add(n));
        lookup.sanityCheck();
        Assert.assertEquals(4, lookup.nextNodes(nodes.get(0)).size());
    }
    @Test
    public void comprehensive_check_nodelength2() {
        for (int k = 3; k <= 4; k++) {
            for (int nodeWidth = 1; nodeWidth <= 4; nodeWidth++) {
                Set<ImmutableKmerNode> existing = new HashSet<>();
                List<KmerNode> nodes = new ArrayList<>();
                for (int startOffset = 0; startOffset <= 6 * nodeWidth; startOffset += nodeWidth + 1) {
                    for (long kmer = 0; kmer < 1 << (2 * k); kmer++) {
                        for (long nextKmer : KmerEncodingHelper.nextStates(k, kmer)) {
                            KmerPathNode kpn = new KmerPathNode(kmer, startOffset, startOffset + nodeWidth - 1, false, nodes.size() + 1);
                            kpn.append(new ImmutableKmerNode(nextKmer, startOffset + 1, startOffset + nodeWidth - 1 + 1, false, nodes.size() + 1));
                            List<ImmutableKmerNode> constituentNodes = constituentNode(kpn);
                            if (Collections.disjoint(existing, constituentNodes)) {
                                nodes.add(kpn);
                                existing.addAll(constituentNodes);
                            }
                        }
                    }
                }
                checkEverything(k, nodes);
            }
        }
    }

    private List<ImmutableKmerNode> constituentNode(KmerPathNode kpn) {
        List<ImmutableKmerNode> list = new ArrayList<>(kpn.length() * kpn.width());
        for (int i = 0; i < kpn.length(); i++) {
            for (int j = kpn.firstStart(); j <= kpn.lastStart(); j++) {
                list.add(new ImmutableKmerNode(kpn.kmer(i), j + i, j + i, false, 0));
            }
        }
        return list;
    }
    @Test
    public void should_handle_loopback_intervals() {
        int k = 4;
        List<ImmutableKmerNode> nodes = ImmutableList.of(
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("CAAA")), 1, 100, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAA")), 2, 101, false, 1),
                new ImmutableKmerNode(KmerEncodingHelper.picardBaseToEncoded(k, B("AAAT")), 3, 102, false, 1));
        final KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
        nodes.stream().forEach(n -> lookup.add(n));
        lookup.sanityCheck();
        Assert.assertEquals(2, lookup.nextNodes(nodes.get(0)).size());
        Assert.assertTrue(lookup.nextNodes(nodes.get(0)).contains(nodes.get(1)));
        Assert.assertTrue(lookup.nextNodes(nodes.get(0)).contains(nodes.get(2)));
        checkEverything(k, nodes);
    }
    public void checkEverything(int k, List<? extends KmerNode> nodes) {
        KmerNodeNonOverlappingLookup<KmerNode> lookup = new KmerNodeNonOverlappingLookup<>(k);
        nodes.stream().forEach(n -> lookup.add(n));
        lookup.sanityCheck();
        for (KmerNode left : nodes) {
            List<KmerNode> nextNodes = lookup.nextNodes(left);
            nextNodes.sort(KmerNodeUtil.ByLastStartEndKmerReferenceWeight);
            List<KmerNode> expectedNodes = nodes.stream()
                    .filter(right ->
                            IntervalUtil.overlapsClosed(left.lastStart() + 1, left.lastEnd() + 1, right.firstStart(), right.firstEnd()) &&
                                    KmerEncodingHelper.isNext(k, left.lastKmer(), right.firstKmer())
                    )
                    .sorted(KmerNodeUtil.ByLastStartEndKmerReferenceWeight)
                    .collect(Collectors.toList());
            Assert.assertEquals(nextNodes, expectedNodes);
            if (nextNodes.size() == 1) {
                KmerNode right = nextNodes.get(0);
                if (left.lastStart() + 1 == right.firstStart() && left.lastEnd() + 1 == right.firstEnd()) {
                    Assert.assertTrue(lookup.getUniqueFullWidthSuccessor(left) == right);
                } else {
                    Assert.assertNull(lookup.getUniqueFullWidthSuccessor(left));
                }
            } else {
                Assert.assertNull(lookup.getUniqueFullWidthSuccessor(left));
            }
        }
        for (KmerNode right : nodes) {
            List<KmerNode> prevNodes = lookup.prevNodes(right);
            prevNodes.sort(KmerNodeUtil.ByLastStartEndKmerReferenceWeight);
            List<KmerNode> expectedNodes = nodes.stream()
                    .filter(left ->
                            IntervalUtil.overlapsClosed(left.lastStart() + 1, left.lastEnd() + 1, right.firstStart(), right.firstEnd()) &&
                                    KmerEncodingHelper.isNext(k, left.lastKmer(), right.firstKmer())
                    )
                    .sorted(KmerNodeUtil.ByLastStartEndKmerReferenceWeight)
                    .collect(Collectors.toList());
            Assert.assertEquals(prevNodes, expectedNodes);
            if (prevNodes.size() == 1) {
                KmerNode left = prevNodes.get(0);
                if (left.lastStart() + 1 == right.firstStart() && left.lastEnd() + 1 == right.firstEnd()) {
                    Assert.assertTrue(left == lookup.getUniqueFullWidthPredecessor(right));
                } else {
                    Assert.assertNull(lookup.getUniqueFullWidthPredecessor(right));
                }
            } else {
                Assert.assertNull(lookup.getUniqueFullWidthPredecessor(right));
            }
        }
        nodes.stream().forEach(n -> lookup.remove(n));
        lookup.sanityCheck();
    }
}
