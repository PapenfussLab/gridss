package au.edu.wehi.idsv.kraken;

import htsjdk.samtools.util.RuntimeIOException;
import joptsimple.internal.Strings;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class KrakenParser implements Closeable, Iterator<KrakenClassification> {
    private final BufferedReader krakenOutput;
    private String nextLine = null;

    public KrakenParser(BufferedReader krakenOutput) {
        this.krakenOutput = krakenOutput;
    }

    @Override
    public void close() throws IOException {
        krakenOutput.close();
    }

    @Override
    public boolean hasNext() {
        if (nextLine == null) {
            try {
                nextLine = krakenOutput.readLine();
            } catch (IOException e) {
                throw new RuntimeIOException(e);
            }
        }
        return !Strings.isNullOrEmpty(nextLine);
    }

    @Override
    public KrakenClassification next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        KrakenClassification kc = new KrakenClassification(nextLine);
        nextLine = null;
        return kc;
    }
}
