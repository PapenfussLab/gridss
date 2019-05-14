package au.edu.wehi.idsv.visualisation;

import java.io.*;
import java.util.Collection;

public class StateTracker implements Closeable {
    private BufferedWriter writer;
    private long lastTime = System.nanoTime();
    private boolean inHeader = true;
    public StateTracker(File file) throws IOException {
        this.writer = new BufferedWriter(new FileWriter(file));
    }

    public void writeHeader(Collection<TrackedState> obj) throws IOException {
        writer.write("nsElapsedTime");
        for (TrackedState o : obj) {
            for (String s : o.trackedNames()) {
                writer.write(",");
                writer.write(s);
            }
        }
        writer.write("\n");
    }

    public void track(Collection<TrackedState> obj) throws IOException {
        for (TrackedState o : obj) {
            for (Object x : o.trackedState()) {
                writer.write(",");
                if (x != null) {
                    writer.write(x.toString());
                }
            }
        }
        writer.write("\n");
    }

    @Override
    public void close() throws IOException {
        writer.close();
    }
}
