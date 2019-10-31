package au.edu.wehi.idsv.util;

import java.util.Iterator;

public class ErrorIterator<T> implements Iterator<T> {

    @Override
    public boolean hasNext() {
        return true;
    }

    @Override
    public T next() {
        throw new RuntimeException("ErrorIterator");
    }
}