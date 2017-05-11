package au.edu.wehi.idsv.visualisation;

import java.io.Closeable;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

import au.edu.wehi.idsv.BreakendDirection;
import htsjdk.samtools.util.Log;

public class AssemblyTelemetry implements Closeable {
	private static final Log log = Log.getInstance(AssemblyTelemetry.class);
	private static HashMap<File, FileWriter> writers = new HashMap<>();
	private static HashMap<File, Integer> counts = new HashMap<>();
	private FileWriter writer;
	private File file;
	private int chunk;
	private BreakendDirection direction;
	private static synchronized FileWriter init(File file) throws IOException {
		if (counts.get(file) != null) {
			counts.put(file, counts.get(file) + 1);
		} else {
			counts.put(file, 1);
			boolean shouldWriteHeader = !file.exists();
			FileWriter fw = new FileWriter(file, true);
			if (shouldWriteHeader) {
				writeHeader(fw);
			}
			writers.put(file, fw);
		}
		return writers.get(file);
	}
	private static synchronized void close(File file) {
		if (counts.get(file) == null) {
			log.debug("Sanity check failure: ", file, " has no open file handles");
			writers.put(file, null);
			return;
		}
		if (counts.get(file) <= 1) {
			counts.remove(file);
			try {
				writers.remove(file).close();
			} catch (IOException e) {
				log.debug(e);
			}
		} else {
			counts.put(file, counts.get(file) - 1);
		}
	}
	private static void writeHeader(FileWriter writer) {
	}
	public AssemblyTelemetry(File telemetryFile, int chunk, BreakendDirection direction) {
		this.file = telemetryFile;
		this.chunk = chunk;
		this.direction = direction;
		try {
			this.writer = init(telemetryFile);
		} catch (IOException e) {
			log.debug(e, String.format("Unable to load chunk %d assembly telemetry for %s", chunk, telemetryFile));
		}
	}
	
	public void loadGraph(int referenceIndex, int start, int end, int nodes, boolean filtered) {
		try {
			writer.write(String.format("%d,%s,load,%d,%d,%d,%d,%b\n", chunk, direction.toChar(), referenceIndex, start, end, nodes, filtered));
		} catch (IOException e) {
		}
	}

	public void flushContigs(int referenceIndex, int flushStart, int flushEnd, int contigsFlushed) {
		try {
			writer.write(String.format("%d,%s,flushContigs,%d,%d,%d,%d,\n", chunk, direction.toChar(), referenceIndex, flushStart, flushEnd, contigsFlushed));
		} catch (IOException e) {
		}
	}

	public void flushReferenceNodes(int referenceIndex, int flushStart, int flushEnd, int readsFlushed) {
		try {
			writer.write(String.format("%d,%s,flushReferenceNodes,%d,%d,%d,%d,\n", chunk, direction.toChar(), referenceIndex, flushStart, flushEnd, readsFlushed));
		} catch (IOException e) {
		}
	}

	public void callContig(int referenceIndex, int start, int end, int nodes, int reads, boolean repeatsSimplified) {
	}

	@Override
	public void close() {
		close(file);
	}
}
