package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.debruijn.positional.KmerSupportNode;
import com.google.common.collect.Iterables;

import java.util.*;

/**
 * Optimised priority queue implementation that efficiently implements only the
 * operations used by SupportNodeIterator
 *
 * Alternate implementation include:
 *  - IntSortedSet lookup backed by Int2ObjectHashMap
 */

public class KmerSupportNodePriorityQueueByFirstStart implements Queue<KmerSupportNode> {
    private final int blockBits;
    private Node head = null;
    private int recordCount = 0;

    public KmerSupportNodePriorityQueueByFirstStart(int blockBits) {
        this.blockBits = blockBits;
    }

    private static class Node {
        public Node next;
        public final int nodeIndex;
        public final ArrayDeque<KmerSupportNode>[] position;
        public int firstOccupiedOffset;
        public int recordCount;
        public Node(int index, int offset, KmerSupportNode n, int blockBits) {
            position = new ArrayDeque[1 << blockBits];
            nodeIndex = index;
            firstOccupiedOffset = offset;
            position[offset] = new ArrayDeque<>();
            position[offset].add(n);
            recordCount = 1;
        }
        public void add(int offset, KmerSupportNode n) {
            if (offset < firstOccupiedOffset) {
                firstOccupiedOffset = offset;
            }
            ArrayDeque<KmerSupportNode> q = position[offset];
            if (q == null) {
                q = new ArrayDeque<>();
                position[offset] = q;
            }
            q.add(n);
            recordCount++;
        }
        private void ensureFirstOccupiedOffsetIsValid() {
            if (recordCount == 0) {
                throw new NoSuchElementException();
            }
            while (position[firstOccupiedOffset] == null) {
                firstOccupiedOffset++;
            }
        }
        public KmerSupportNode popFirst() {
            ensureFirstOccupiedOffsetIsValid();
            recordCount--;
            KmerSupportNode n = position[firstOccupiedOffset].removeFirst();
            if (position[firstOccupiedOffset].size() == 0) {
                position[firstOccupiedOffset] = null;
            }
            return n;
        }
        public KmerSupportNode peekFirst() {
            ensureFirstOccupiedOffsetIsValid();
            return position[firstOccupiedOffset].peekFirst();
        }
    }

    @Override
    public int size() {
        return recordCount;
    }

    @Override
    public boolean isEmpty() {
        return recordCount == 0;
    }

    @Override
    public boolean contains(Object o) {
        throw new UnsupportedOperationException();
    }

    @Override
    public Iterator<KmerSupportNode> iterator() {
        throw new UnsupportedOperationException();
    }

    @Override
    public Object[] toArray() {
        throw new UnsupportedOperationException();
    }

    @Override
    public <T> T[] toArray(T[] a) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean add(KmerSupportNode kmerSupportNode) {
        recordCount++;
        int position = kmerSupportNode.firstStart();
        int nodeIndex = position >> blockBits;
        int nodeOffset = position - (nodeIndex << blockBits);
        if (head == null) {
            head = new Node(nodeIndex, nodeOffset, kmerSupportNode, blockBits);
        } else if (head.nodeIndex > nodeIndex) {
            // prepend this record before the current head
            Node oldHead = head;
            head = new Node(nodeIndex, nodeOffset, kmerSupportNode, blockBits);
            head.next = oldHead;
        } else {
            Node node = head;
            while (node.next != null && node.next.nodeIndex <= nodeIndex) {
                node = node.next;
            }
            if (node.nodeIndex == nodeIndex) {
                // Found our node
                node.add(nodeOffset, kmerSupportNode);
            } else {
                // node after the current one
                Node oldNext = node.next;
                node.next = new Node(nodeIndex, nodeOffset, kmerSupportNode, blockBits);
                node.next.next = oldNext;
            }
        }
        //assert(sanityCheck());
        return true;
    }

    @Override
    public boolean remove(Object o) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean addAll(Collection<? extends KmerSupportNode> c) {
        for (KmerSupportNode ksn : c) {
            add(ksn);
        }
        return true;
    }

    @Override
    public boolean removeAll(Collection<?> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean retainAll(Collection<?> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void clear() {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean offer(KmerSupportNode kmerSupportNode) {
        throw new UnsupportedOperationException();
    }

    @Override
    public KmerSupportNode remove() {
        throw new UnsupportedOperationException();
    }

    @Override
    public KmerSupportNode poll() {
        KmerSupportNode n = head.popFirst();
        if (head.recordCount == 0) {
            head = head.next;
        }
        recordCount--;
        //assert(sanityCheck());
        return n;
    }

    @Override
    public KmerSupportNode element() {
        throw new UnsupportedOperationException();
    }

    @Override
    public KmerSupportNode peek() {
        //assert(sanityCheck());
        return head.peekFirst();
    }
    public boolean sanityCheck() {
        Node n = head;
        while (n != null && n.next != null) {
            assert(n.nodeIndex < n.next.nodeIndex);
        }
        n = head;
        int totalCount = 0;
        while (n != null) {
            int nodeCount = 0;
            for (int i = 0; i < n.position.length; i++) {
                ArrayDeque<KmerSupportNode> q = n.position[i];
                if (q != null) {
                    nodeCount += q.size();
                    final Node nn = n;
                    final int ii = i;
                    Iterables.all(q, x -> x.firstStart() == nn.nodeIndex * (1 << blockBits) + ii);
                    Set<KmerSupportNode> s = new HashSet<>(q);
                    assert(q.size() == s.size());
                }
            }
            assert(nodeCount == n.recordCount);
            n = n.next;
        }
        assert(totalCount == recordCount);
        return true;
    }
}
