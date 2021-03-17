package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.positional.ImmutableKmerNode;
import au.edu.wehi.idsv.debruijn.positional.KmerNode;
import au.edu.wehi.idsv.util.IntervalUtil;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

public class KmerNodeByFirstKmerIntervalLookupTest extends TestHelper {
    private static final int k = 4;
    private KmerNode kn(String kmer, int start, int end) {
        return kn(KmerEncodingHelper.picardBaseToEncoded(k, kmer.getBytes()), start, end);
    }
    private KmerNode kn(long kmer, int start, int end) {
        return new ImmutableKmerNode(kmer, start, end, false, 1);
    }
    @Test
    public void get_should_match_overlap_logic() {
        for (List<KmerNode> kns : ImmutableList.of(
                ImmutableList.of(
                        kn(0, 1, 2),
                        kn(0, 4, 5),
                        kn(0, 6, 6),
                        kn(0, 7, 7),
                        kn(0, 9, 9)),
                ImmutableList.of(
                        kn(0, 1, 3),
                        kn(0, 7, 7)),
                ImmutableList.of(
                        kn(0, 1, 2)),
                ImmutableList.of(
                        kn(0, 1, 2),
                        kn(1, 1, 2))
        )) {
            KmerNodeByLastKmerIntervalLookup<KmerNode> lookup = new KmerNodeByLastKmerIntervalLookup<>();
            kns.stream().forEach(n -> lookup.add(n));
            validate_against_direct_comparison(lookup, kns);
        }
    }

    private void validate_against_direct_comparison(KmerNodeByLastKmerIntervalLookup<KmerNode> lookup, List<KmerNode> kns) {
        for (long kmer : new long[] { 0, 1, 2}) {
            for (int i = -1; i < 11; i++) {
                for (int j = i; j < 12; j++) {
                    int start = i;
                    int end = j;
                    List<KmerNode> expected = kns.stream()
                            .filter(n -> n.firstKmer() == kmer && IntervalUtil.overlapsClosed(start, end, n.firstStart(), n.firstEnd()))
                            .collect(Collectors.toList());
                    Assert.assertEquals(expected, lookup.getOverlapping(kmer, i, j));
                    KmerNode exactMatch = kns.stream()
                            .filter(n -> n.firstKmer() == kmer && n.firstStart() == start && n.firstEnd() == end)
                            .findFirst().orElse(null);
                    Assert.assertEquals(exactMatch, lookup.get(kmer, i, j));
                }
            }
        }
    }

    @Test
    public void remove_should_match_overlap_logic() {
        for (List<KmerNode> full : ImmutableList.of(
                ImmutableList.of(
                        kn(0, 1, 2),
                        kn(0, 4, 5),
                        kn(0, 6, 6),
                        kn(0, 7, 7),
                        kn(0, 9, 9)),
                ImmutableList.of(
                        kn(0, 1, 3),
                        kn(0, 7, 7)),
                ImmutableList.of(
                        kn(0, 1, 2)),
                ImmutableList.of(
                        kn(0, 1, 2),
                        kn(1, 1, 2))
        )) {
            for (int offset = 0; offset < full.size(); offset++) {
                ArrayList<KmerNode> kns = Lists.newArrayList(full);
                KmerNodeByLastKmerIntervalLookup<KmerNode> lookup = new KmerNodeByLastKmerIntervalLookup<>();
                full.stream().forEach(n -> lookup.add(n));
                // remove each element
                KmerNode removed = kns.remove(offset);
                lookup.remove(removed);
                validate_against_direct_comparison(lookup, kns);
                for (int offset2 = 0; offset2 < kns.size(); offset2++) {
                    KmerNodeByLastKmerIntervalLookup lookup2 = new KmerNodeByLastKmerIntervalLookup<>();
                    full.stream().forEach(n -> lookup2.add(n));
                    ArrayList<KmerNode> kns2 = Lists.newArrayList(kns);
                    KmerNode removed2 = kns2.remove(offset2);
                    lookup2.remove(removed);
                    lookup2.remove(removed2);
                    validate_against_direct_comparison(lookup2, kns2);
                }
            }
        }
    }
}