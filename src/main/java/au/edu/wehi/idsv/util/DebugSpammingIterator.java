package au.edu.wehi.idsv.util;

import htsjdk.samtools.util.RuntimeIOException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;

public class DebugSpammingIterator<T> implements Iterator<T> {
    private static final int FLUSH_INTERVAL = 1024;
    private static final HashMap<String, PrintWriter> out = new HashMap<>();
    private final Iterator<T> underlying;
    private final long maxRecords;
    private final String id;
    private long records = 0;
    public DebugSpammingIterator(Iterator<T> underlying, String id) {
        this(underlying, id, 1024 * 1024 * 1024);
    }
    public DebugSpammingIterator(Iterator<T> underlying, String id, long maxRecords) {
        this.underlying = underlying;
        this.id = id;
        this.maxRecords = maxRecords;
    }
    @Override
    public boolean hasNext() {
        return underlying.hasNext();
    }

    @Override
    public T next() {
        T r = underlying.next();
        String threadName = Thread.currentThread().getName();
        records++;
        if (records < maxRecords) {
            if (!out.containsKey(threadName)) {
                int i = 1;
                while (new File("DebugSpammingIterator" + i + "." + threadName + ".log").exists()) {
                    i++;
                }
                try {
                    out.put(threadName, new PrintWriter(new File("DebugSpammingIterator" + i + "." + threadName + ".log")));
                } catch (FileNotFoundException e) {
                    throw new RuntimeIOException(e);
                }
            }
            String msg = String.format("%s\t%s\t%s\n", Thread.currentThread().getName(), id, r);
            out.get(threadName).write(msg);
        } else if (records == maxRecords || records % FLUSH_INTERVAL == 0) {
            synchronized (out) {
                out.get(threadName).flush();
            }
        }
        return r;
    }
}
