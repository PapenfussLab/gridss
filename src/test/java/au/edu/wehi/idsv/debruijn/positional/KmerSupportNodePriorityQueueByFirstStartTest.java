package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.SingleReadEvidence;
import au.edu.wehi.idsv.TestHelper;
import org.junit.Ignore;
import org.junit.Test;

import java.util.PriorityQueue;
import java.util.Random;

import static org.junit.Assert.*;

public class KmerSupportNodePriorityQueueByFirstStartTest extends TestHelper {
    @Test
    public void should_mimic_priority_queue() {
        PriorityQueue<KmerSupportNode> pq = new PriorityQueue<>(1024, KmerNodeUtil.ByFirstStart);
        KmerSupportNodePriorityQueueByFirstStart optimisepq = new KmerSupportNodePriorityQueueByFirstStart(4);
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
        while (!pq.isEmpty()) {
            assertFalse(optimisepq.isEmpty());
            assertEquals(pq.peek(), optimisepq.peek());
            assertEquals(pq.poll(), optimisepq.poll());
            assertEquals(pq.size(), optimisepq.size());
        }
        for (int i = 9000; i >= 1; i--) {
            KmerSupportNode ksn = ksnlist[i];
            pq.add(ksn);
            optimisepq.add(ksn);
            assertEquals(pq.size(), optimisepq.size());
            assertEquals(pq.peek(), optimisepq.peek());
        }
        while (!pq.isEmpty()) {
            assertFalse(optimisepq.isEmpty());
            assertEquals(pq.peek(), optimisepq.peek());
            assertEquals(pq.poll(), optimisepq.poll());
            assertEquals(pq.size(), optimisepq.size());
        }
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
        for (int i = 0; i < 100000; i++) {
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
        while (!pq.isEmpty()) {
            assertFalse(optimisepq.isEmpty());
            assertEquals(pq.peek(), optimisepq.peek());
            assertEquals(pq.poll(), optimisepq.poll());
            assertEquals(pq.size(), optimisepq.size());
        }
    }
    @Test
    @Ignore("Class not yet in use")
    public void should_handle_negative_positions() {
        PriorityQueue<KmerSupportNode> pq = new PriorityQueue<>(1024, KmerNodeUtil.ByFirstStart);
        KmerSupportNodePriorityQueueByFirstStart optimisepq = new KmerSupportNodePriorityQueueByFirstStart(4);
        MockSAMEvidenceSource ses = SES();
        KmerSupportNode[] ksnlist = new KmerSupportNode[9002];
        for (int i = 1; i < 100; i++) {
            KmerEvidence evidence = KmerEvidence.create(1, NRRP(ses, DP(0, i, "10M", false, 1, i, "10M", false)));
            ksnlist[i] = evidence.node(0);
            KmerSupportNode ksn = ksnlist[i];
            pq.add(ksn);
            optimisepq.add(ksn);
            assertEquals(pq.size(), optimisepq.size());
            assertEquals(pq.peek(), optimisepq.peek());
        }
        while (!pq.isEmpty()) {
            assertFalse(optimisepq.isEmpty());
            assertEquals(pq.peek(), optimisepq.peek());
            assertEquals(pq.poll(), optimisepq.poll());
            assertEquals(pq.size(), optimisepq.size());
        }
    }
}