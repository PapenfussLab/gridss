package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.debruijn.positional.KmerPathNode;
import au.edu.wehi.idsv.debruijn.positional.KmerPathSubnode;
import au.edu.wehi.idsv.debruijn.positional.TraversalNode;
import com.google.common.collect.Lists;
import org.junit.Test;

import java.util.ArrayList;
import java.util.SortedSet;
import java.util.TreeSet;

import static org.junit.Assert.assertEquals;

public class TraversalNodeByScoreDescPathFirstIdentityTest {
    @Test
    public void should_match_sortedset() {
        ArrayList<TraversalNode> list = new ArrayList<>();
        for(int score = 0; score < 4; score++) {
            for (int pathStart = 0; pathStart < 4; pathStart++) {
                for (int pathEnd = pathStart; pathEnd <= pathStart + 1; pathEnd++) {
                    KmerPathNode kpn = new KmerPathNode(list.size(), pathStart, pathEnd, false, score);
                    KmerPathSubnode kps = new KmerPathSubnode(kpn, pathStart, pathEnd);
                    list.add(new TraversalNode(kps, score));
                }
            }
        }
        SortedSet<TraversalNode> ss = new TreeSet<>(TraversalNode.ByScoreDescPathFirstEndSubnode);
        TraversalNodeByScoreDescPathFirstIdentity opt = new TraversalNodeByScoreDescPathFirstIdentity();
        for (TraversalNode tn : list) {
            assertEquals(ss.add(tn), opt.add(tn));
            assertEquals(ss.size(), opt.size());
            assertEquals(ss.add(tn), opt.add(tn));
            assertEquals(ss.size(), opt.size());
            assertEquals(ss.remove(tn), opt.remove(tn));
            assertEquals(ss.size(), opt.size());
            assertEquals(ss.add(tn), opt.add(tn));
            assertEquals(ss.size(), opt.size());
            assertEquals(ss.first(), opt.first());
        }
        assertEquals(ss.removeAll(list), opt.removeAll(list));
        assertEquals(ss.size(), opt.size());
        for (TraversalNode tn : Lists.reverse(list)) {
            assertEquals(ss.add(tn), opt.add(tn));
            assertEquals(ss.first(), opt.first());
        }
    }
}