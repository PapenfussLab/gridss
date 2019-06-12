package au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures;

import au.edu.wehi.idsv.debruijn.positional.KmerEvidence;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

public class KmerEvidenceSet implements Set<KmerEvidence> {
    private HashMap<String, KmerEvidence> lookup = new HashMap<>();
    private String lookupKey(KmerEvidence e) {
        return e.evidence().getEvidenceID();
    }
    @Override
    public int size() {
        return lookup.size();
    }

    @Override
    public boolean isEmpty() {
        return lookup.isEmpty();
    }

    @Override
    public boolean contains(Object o) {
        if (o instanceof KmerEvidence) {
            return lookup.containsKey(lookupKey((KmerEvidence)o));
        }
        return false;
    }

    @Override
    public Iterator<KmerEvidence> iterator() {
        return lookup.values().iterator();
    }

    @Override
    public Object[] toArray() {
        return lookup.values().toArray();
    }

    @Override
    public <T> T[] toArray(T[] a) {
        return lookup.values().toArray(a);
    }

    @Override
    public boolean add(KmerEvidence kmerEvidence) {
        return lookup.put(lookupKey(kmerEvidence), kmerEvidence) == null;
    }

    @Override
    public boolean remove(Object o) {
        if (o instanceof KmerEvidence) {
            return lookup.remove(lookupKey((KmerEvidence)o)) != null;
        }
        return false;
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        for (Object o : c) {
            if (!contains(o)) return false;
        }
        return true;
    }

    @Override
    public boolean addAll(Collection<? extends KmerEvidence> c) {
        boolean changed = false;
        for (KmerEvidence e : c) {
            changed |= add(e);
        }
        return changed;
    }

    @Override
    public boolean retainAll(Collection<?> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean removeAll(Collection<?> c) {
        boolean changed = false;
        for (Object o : c) {
            changed |= remove(o);
        }
        return changed;
    }

    @Override
    public void clear() {
        lookup.clear();
    }
}
