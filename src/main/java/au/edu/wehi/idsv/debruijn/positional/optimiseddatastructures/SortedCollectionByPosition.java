package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import java.util.NoSuchElementException;

/**
 * Collection data structure keyed by position optimised
 * for windowed linear traversal through the genome
 *
 * Alternate implementation include:
 *  - IntSortedSet lookup backed by Int2ObjectHashMap
 */
public abstract class SortedCollectionByPosition<T, TColl> {
    private final int blockBits;
    private Node<TColl> head = null;
    private int recordCount = 0;

    public SortedCollectionByPosition(int blockBits) {
        this.blockBits = blockBits;
    }

    private static class Node<TColl> {
        public Node next;
        public final int nodeIndex;
        public final TColl[] position;
        public int firstOccupiedOffset;
        public int recordCount;
        public Node(int index, int blockBits) {
            position = (TColl[])new Object[1 << blockBits];
            nodeIndex = index;
            recordCount = 1;
        }
    }
    private void ensureNodeFirstOccupiedOffsetIsValid(Node n) {
        if (n.recordCount == 0) {
            throw new NoSuchElementException();
        }
        while (n.position[n.firstOccupiedOffset] == null) {
            n.firstOccupiedOffset++;
        }
    }
    private boolean addToNode(Node<TColl> n, int offset, T obj) {
        this.recordCount++;
        n.recordCount++;
        if (offset < n.firstOccupiedOffset) {
            n.firstOccupiedOffset = offset;
        }
        if (n.position[offset] == null) {
            n.position[offset] = createAtPosition();
        }
        return addAtPosition(n.position[offset], obj);
    }
    private Node<TColl> createNode(int index, int offset, T obj) {
        Node<TColl> n = new Node<>(index, blockBits);
        n.firstOccupiedOffset = offset;
        n.position[offset] = createAtPosition();
        addAtPosition(n.position[offset], obj);
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
    protected abstract boolean dataIsEmpty(TColl coll);

    public int size() {
        return recordCount;
    }

    public boolean isEmpty() {
        return recordCount == 0;
    }

    public boolean add(T x) {
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
        recordCount++;
        return true;
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
        if (dataIsEmpty(node.position[nodeOffset])) {
            node.position[nodeOffset] = null;
        }
        node.recordCount--;
        if (node.recordCount == 0) {
            if (prevNode == null) {
                head = node.next;
            } else {
                prevNode.next = node.next;
            }
        }
        return result;
    }

    public void clear() {
        recordCount = 0;
        head = null;
    }

    protected T pollFirstEntry() {
        if (head == null) throw new NoSuchElementException();
        ensureNodeFirstOccupiedOffsetIsValid(head);
        T result = popAtPosition(head.position[head.firstOccupiedOffset]);
        if (dataIsEmpty(head.position[head.firstOccupiedOffset])) {
            head.position[head.firstOccupiedOffset] = null;
        }
        head.recordCount--;
        if (head.recordCount == 0) {
            head = head.next;
        }
        return result;
    }

    protected T peekFirstEntry() {
        ensureNodeFirstOccupiedOffsetIsValid(head);
        return peekAtPosition(head.position[head.firstOccupiedOffset]);
    }
}
