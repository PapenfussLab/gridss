package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.stream.Stream;

/**
 * Collection data structure keyed by position optimised
 * for windowed linear traversal through the genome
 *
 * Alternate implementation include:
 *  - IntSortedSet lookup backed by Int2ObjectHashMap
 */
public abstract class SortedByPosition<T, TColl> {
    private final int blockBits;
    private Node<TColl> head = null;
    private int recordCount = 0;

    public SortedByPosition(int blockBits) {
        this.blockBits = blockBits;
    }

    private static class Node<TColl> {
        public final int nodeIndex;
        public final TColl[] position;
        public Node next = null;
        public int firstOccupiedOffset = 0;
        public int recordCount = 0;
        public Node(int index, int blockBits) {
            position = (TColl[])new Object[1 << blockBits];
            nodeIndex = index;
        }
    }
    private void ensureHeadFirstOccupiedOffsetIsValid() {
        while (head != null && head.recordCount == 0) {
            head = head.next;
        }
        if (head == null) {
            throw new NoSuchElementException();
        }
        while(head.firstOccupiedOffset < head.position.length) {
            TColl t = head.position[head.firstOccupiedOffset];
            if (t != null) {
                if (positionIsEmpty(t)) {
                    head.position[head.firstOccupiedOffset] = null;
                } else {
                    return;
                }
            }
            head.firstOccupiedOffset++;
        }
        throw new IllegalStateException("Sanity check failure: overrun array bounds");
    }
    private boolean addToNode(Node<TColl> n, int offset, T obj) {
        if (offset < n.firstOccupiedOffset) {
            n.firstOccupiedOffset = offset;
        }
        if (n.position[offset] == null) {
            TColl coll = createAtPosition();
            assert(coll != null);
            n.position[offset] = coll;
        }
        boolean added = addAtPosition(n.position[offset], obj);
        if (added) {
            this.recordCount++;
            n.recordCount++;
        }
        return added;
    }
    private Node<TColl> createNode(int index, int offset, T obj) {
        Node<TColl> n = new Node<>(index, blockBits);
        n.firstOccupiedOffset = offset;
        addToNode(n, offset, obj);
        return n;
    }

    /**
     * Gets the genomic positional key
     */
    protected abstract int getPosition(T obj);

    /**
     * Gets the record at the first genomic position
     */
    protected abstract T peekAtPosition(TColl coll);

    /**
     * Removes the head element for the given genomic position
     */
    protected abstract T popAtPosition(TColl coll);

    /**
     * Creates a new collection for a genomic position
     */
    protected abstract TColl createAtPosition();

    /**
     * Adds the given element to the given genomic position collection
     */
    protected abstract boolean addAtPosition(TColl existing, T toAdd);

    /**
     * Removes the given element at the collection for a given genomic position
     * @return returns a tuple indicating whether the collection is now empty
     * and whether the given element was removed, in that order
     */
    protected abstract boolean removeAtPosition(TColl coll, T obj);

    /**
     * @return true if the given data structure is empty, false otherwise.
     */
    protected abstract boolean positionIsEmpty(TColl coll);

    protected abstract boolean containsAtPosition(TColl coll, T obj);

    public int size() {
        return recordCount;
    }

    public boolean isEmpty() {
        return recordCount == 0;
    }

    public boolean add(T x) {
        int preCount = recordCount;
        int position = getPosition(x);
        int nodeIndex = position >> blockBits;
        int nodeOffset = position - (nodeIndex << blockBits);
        if (head == null) {
            head = createNode(nodeIndex, nodeOffset, x);
        } else if (head.nodeIndex > nodeIndex) {
            // prepend this record before the current head
            Node<TColl> oldHead = head;
            head = createNode(nodeIndex, nodeOffset, x);
            head.next = oldHead;
        } else {
            Node<TColl> node = head;
            while (node.next != null && node.next.nodeIndex <= nodeIndex) {
                node = node.next;
            }
            if (node.nodeIndex == nodeIndex) {
                // found our node
                addToNode(node, nodeOffset, x);
            } else {
                // add new node after the current one
                Node oldNext = node.next;
                node.next = createNode(nodeIndex, nodeOffset, x);
                node.next.next = oldNext;
            }
        }
        assert(sanityCheck());
        return recordCount != preCount;
    }

    public boolean remove(Object x) {
        T t = (T)x;
        int position = getPosition(t);
        int nodeIndex = position >> blockBits;
        int nodeOffset = position - (nodeIndex << blockBits);
        Node<TColl> prevNode = null;
        Node<TColl> node = head;
        while (node != null && node.nodeIndex < nodeIndex) {
            prevNode = node;
            node = node.next;
        }
        if (node == null || node.nodeIndex > nodeIndex) {
            return false;
        }
        if (node.position[nodeOffset] == null) {
            return false;
        }
        boolean result = removeAtPosition(node.position[nodeOffset], t);
        if (result) {
            if (positionIsEmpty(node.position[nodeOffset])) {
                node.position[nodeOffset] = null;
            }
            recordCount--;
            node.recordCount--;
            if (node.recordCount == 0) {
                if (prevNode == null) {
                    head = node.next;
                } else {
                    prevNode.next = node.next;
                }
            }
        }
        assert(sanityCheck());
        return result;
    }

    public void clear() {
        recordCount = 0;
        head = null;
        assert(sanityCheck());
    }

    protected T removeFirstEntry(boolean throwIfEmpty) {
        if (isEmpty()) {
            if (throwIfEmpty) {
                throw new NoSuchElementException();
            } else {
                return null;
            }
        }
        ensureHeadFirstOccupiedOffsetIsValid();
        T result = popAtPosition(head.position[head.firstOccupiedOffset]);
        if (positionIsEmpty(head.position[head.firstOccupiedOffset])) {
            head.position[head.firstOccupiedOffset] = null;
        }
        if (result != null) {
            recordCount--;
            head.recordCount--;
            if (head.recordCount == 0) {
                head = head.next;
            }
        }
        assert(sanityCheck());
        return result;
    }

    protected T peekFirstEntry(boolean throwIfEmpty) {
        if (isEmpty()) {
            if (throwIfEmpty) {
                throw new NoSuchElementException();
            } else {
                return null;
            }
        }
        ensureHeadFirstOccupiedOffsetIsValid();
        return peekAtPosition(head.position[head.firstOccupiedOffset]);
    }
    public boolean contains(Object x) {
        T t = (T)x;
        int position = getPosition(t);
        int nodeIndex = position >> blockBits;
        int nodeOffset = position - (nodeIndex << blockBits);
        Node<TColl> node = head;
        while (node != null && node.nodeIndex < nodeIndex) {
            node = node.next;
        }
        if (node != null && node.nodeIndex == nodeIndex) {
            if (node.position[nodeOffset] == null) return false;
            return containsAtPosition(node.position[nodeOffset], t);
        }
        return false;
    }

    public T remove() {
        return removeFirstEntry(true);
    }

    public T poll() {
        return removeFirstEntry(false);
    }

    public T element() {
        return peekFirstEntry(true);
    }

    public T peek() {
        return peekFirstEntry(false);
    }

    public T first() {
        return peekFirstEntry(true);
    }

    public Object[] toArray()  {
        throw new UnsupportedOperationException();
    }

    public <T> T[] toArray(T[] a)  {
        throw new UnsupportedOperationException();
    }

    public boolean containsAll(Collection<?> c) {
        for (Object o : c) {
            if (!contains(o)) {
                return false;
            }
        }
        return true;
    }

    public boolean addAll(Collection<? extends T> c) {
        boolean changed = false;
        for (T o : c) {
            changed |= add(o);
        }
        return changed;
    }

    public boolean removeAll(Collection<?> c) {
        boolean changed = false;
        for (Object o : c) {
            changed |= remove(o);
        }
        assert(sanityCheck());
        return changed;
    }

    public boolean retainAll(Collection<?> c)  {
        throw new UnsupportedOperationException();
    }
    public boolean sanityCheck() {
        Node n = head;
        int count = 0;
        while (n != null) {
            if (n.next != null) {
                assert (n.nodeIndex < n.next.nodeIndex);
            }
            assert(n.recordCount > 0);
            count += n.recordCount;
            assert(n.recordCount == 0 || Stream.of(n.position).anyMatch(x -> x != null && !positionIsEmpty((TColl)x)));
            n = n.next;
        }
        assert(recordCount == count);
        return true;
    }

    protected List<TColl> getCollectionsInOrder() {
        List<TColl> list = new ArrayList<>();
        for (Node n = head; n != null; n = n.next) {
            for (int i = 0; i < n.position.length; i++) {
                if (n.position[i] != null) {
                    list.add((TColl)n.position[i]);
                }
            }
        }
        return list;
    }
}
