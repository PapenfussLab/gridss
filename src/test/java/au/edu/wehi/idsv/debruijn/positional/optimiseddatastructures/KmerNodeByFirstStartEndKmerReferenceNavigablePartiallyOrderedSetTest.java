package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.positional.KmerNodeUtil;
import au.edu.wehi.idsv.debruijn.positional.KmerPathNode;
import org.junit.Test;

import java.util.NavigableSet;
import java.util.Random;
import java.util.TreeSet;

import static org.junit.Assert.assertEquals;

public class KmerNodeByFirstStartEndKmerReferenceNavigablePartiallyOrderedSetTest extends TestHelper {
    @Test
    public void should_match_navigable_set() {
        int k = 4;
        KmerPathNode[] list = new KmerPathNode[] {
            KPN(k, "GTAC", 1, 10, false),
            KPN(k, "GTAC", 1, 10, true),
            KPN(k, "GTAC", 1, 9, true),
            KPN(k, "GTAC", 0, 10, true),
            KPN(k, "GTAC", 1, 11, true),
            KPN(k, "TTAC", 1, 10, true),
        };
        Random r = new Random(0);
        NavigableSet<KmerPathNode> ns = new TreeSet<>(KmerNodeUtil.ByFirstStartEndKmerReference);
        KmerNodeByFirstStartEndKmerReferenceNavigablePartiallyOrderedSet set = new KmerNodeByFirstStartEndKmerReferenceNavigablePartiallyOrderedSet(4);
        for (int i = 0 ; i < 4096; i++) {
            KmerPathNode kpn = list[r.nextInt(list.length)];
            assertEquals(ns.contains(kpn), set.contains(kpn));
            if (r.nextInt(5) < 2) {
                assertEquals(ns.remove(kpn), set.remove(kpn));
            } else {
                assertEquals(ns.add(kpn), set.add(kpn));
            }
            assertEquals(ns.size(), set.size());
            assertEquals(ns.contains(kpn), set.contains(kpn));
            if (!ns.isEmpty()) {
                // TODO: do we require a stable ordering?
                //assertEquals(ns.first(), set.first());
            }
        }
    }
}