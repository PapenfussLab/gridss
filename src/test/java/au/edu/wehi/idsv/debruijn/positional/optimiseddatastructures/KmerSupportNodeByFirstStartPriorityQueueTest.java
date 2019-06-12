package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.SingleReadEvidence;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.positional.KmerEvidence;
import au.edu.wehi.idsv.debruijn.positional.KmerNodeUtil;
import au.edu.wehi.idsv.debruijn.positional.KmerSupportNode;
import org.junit.Test;

import java.util.PriorityQueue;
import java.util.Random;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

public class KmerSupportNodeByFirstStartPriorityQueueTest extends TestHelper {
    @Test
    public void should_mimic_priority_queue() {
        PriorityQueue<KmerSupportNode> pq = new PriorityQueue<>(1024, KmerNodeUtil.ByFirstStart);
        KmerSupportNodeByFirstStartPriorityQueue optimisepq = new KmerSupportNodeByFirstStartPriorityQueue(4);
        MockSAMEvidenceSource ses = SES();
        KmerSupportNode[] ksnlist = new KmerSupportNode[9002];
        for (int i = 0; i < 9001; i++) {
            KmerEvidence evidence = KmerEvidence.create(1, SingleReadEvidence.createEvidence(ses, 0, Read(0, i, "1M1S")).get(0));
            ksnlist[i] = evidence.node(0);
        }
        for (int i = 1; i < 9000; i++) {
            KmerSupportNode ksn = ksnlist[i];
            pq.add(ksn);
            optimisepq.add(ksn);
            assertEquals(pq.size(), optimisepq.size());
            assertEquals(pq.peek(), optimisepq.peek());
        }
        flush(pq, optimisepq);
        for (int i = 9000; i >= 1; i--) {
            KmerSupportNode ksn = ksnlist[i];
            pq.add(ksn);
            optimisepq.add(ksn);
            assertEquals(pq.size(), optimisepq.size());
            assertEquals(pq.peek(), optimisepq.peek());
        }
        flush(pq, optimisepq);
        for (int i = 1; i < 9000; i++) {
            KmerSupportNode ksn = ksnlist[i];
            pq.add(ksn);
            optimisepq.add(ksn);
            assertEquals(pq.size(), optimisepq.size());
            assertEquals(pq.peek(), optimisepq.peek());
            assertEquals(pq.poll(), optimisepq.poll());
            assertEquals(pq.size(), optimisepq.size());
        }
        Random rng = new Random(0);
        for (int i = 0; i < 10000; i++) {
            KmerSupportNode ksn = ksnlist[1 + rng.nextInt(8999)];
            pq.add(ksn);
            optimisepq.add(ksn);
            assertEquals(pq.size(), optimisepq.size());
            assertEquals(pq.peek(), optimisepq.peek());
            if (rng.nextBoolean()) {
                assertEquals(pq.poll(), optimisepq.poll());
                assertEquals(pq.size(), optimisepq.size());
            }
        }
        flush(pq, optimisepq);
    }
    @Test
    public void should_handle_negative_positions() {
        PriorityQueue<KmerSupportNode> pq = new PriorityQueue<>(1024, KmerNodeUtil.ByFirstStart);
        KmerSupportNodeByFirstStartPriorityQueue optimisepq = new KmerSupportNodeByFirstStartPriorityQueue(3);
        MockSAMEvidenceSource ses = SES();
        KmerSupportNode[] ksnlist = new KmerSupportNode[500];
        for (int i = 0; i < 500; i++) {
            ksnlist[i] = KmerEvidence.create(1, NRRP(ses, DP(0, i, "10M", false, 1, i+1, "10M", false))).node(0);
        }
        Random rng = new Random(0);
        for (int i = 0; i < 10000; i++) {
            KmerSupportNode ksn = ksnlist[rng.nextInt(ksnlist.length)];
            pq.add(ksn);
            optimisepq.add(ksn);
            assertEquals(pq.size(), optimisepq.size());
            assertEquals(pq.peek(), optimisepq.peek());
            if (rng.nextBoolean()) {
                assertEquals(pq.poll(), optimisepq.poll());
                assertEquals(pq.size(), optimisepq.size());
            }
        }
        flush(pq, optimisepq);
    }

    private void flush(PriorityQueue<KmerSupportNode> pq, KmerSupportNodeByFirstStartPriorityQueue optimisepq) {
        while (!pq.isEmpty()) {
            assertFalse(optimisepq.isEmpty());
            assertEquals(pq.peek(), optimisepq.peek());
            assertEquals(pq.poll(), optimisepq.poll());
            assertEquals(pq.size(), optimisepq.size());
        }
    }
}
