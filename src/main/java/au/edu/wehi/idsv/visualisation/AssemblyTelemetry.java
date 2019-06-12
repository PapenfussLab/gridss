package au.edu.wehi.idsv.visualisation;

import au.edu.wehi.idsv.BreakendDirection;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Log;

import java.io.Closeable;
import java.io.File;
import java.io.FileWriter;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

public class AssemblyTelemetry implements Closeable {
	private static final Log log = Log.getInstance(AssemblyTelemetry.class);
	private final File file;
	private final SAMSequenceDictionary dict;
	private BlockingQueue<String> queue;
	public AssemblyTelemetry(File telemetryFile, SAMSequenceDictionary dict) {
		this.file = telemetryFile;
		this.queue = new ArrayBlockingQueue<>(4096);
		this.dict = dict;
		Thread thread = new Thread(new WriterRunnable(), "AT:" + file.getName());
		thread.setDaemon(true);
		thread.start();
	}
	private static void writeHeader(FileWriter writer) {
	}
	public AssemblyChunkTelemetry getTelemetry(int chunkNumber, BreakendDirection direction) {
		return new AssemblyChunkTelemetry(chunkNumber, direction);
	}
	public class AssemblyChunkTelemetry {
		private int chunk;
		private BreakendDirection direction;
		private AssemblyChunkTelemetry(int chunk, BreakendDirection direction) {
			this.chunk = chunk;
			this.direction = direction;
		}
		public void loadGraph(int referenceIndex, int start, int end, int nodes, boolean filtered, long nsSinceLast) {
			String str = String.format("%d,%s,load,%s,%d,%d,%d,%b,%d\n", chunk, direction.toChar(), dict.getSequence(referenceIndex).getSequenceName(), start, end, nodes, filtered, nsSinceLast / 1000);
			put(str);
		}

		public void flushContigs(int referenceIndex, int flushStart, int flushEnd, int contigsFlushed, long nsSinceLast) {
			String str = String.format("%d,%s,flushContigs,%s,%d,%d,%d,,%d\n", chunk, direction.toChar(), dict.getSequence(referenceIndex).getSequenceName(), flushStart, flushEnd, contigsFlushed, nsSinceLast / 1000);
			put(str);
		}

		public void flushReferenceNodes(int referenceIndex, int flushStart, int flushEnd, int readsFlushed, long nsSinceLast) {
			String str = String.format("%d,%s,flushReferenceNodes,%s,%d,%d,%d,,%d\n", chunk, direction.toChar(), dict.getSequence(referenceIndex).getSequenceName(), flushStart, flushEnd, readsFlushed, nsSinceLast / 1000);
			put(str);
		}
		public void callContig(int referenceIndex, int start, int end, int nodes, int reads, boolean repeatsSimplified) {
		}
	}
	private void put(String str) {
		try {
			if (queue != null) {
				queue.put(str);
			}
		} catch (InterruptedException e) {
			queue = null;
		}
	}
	@Override
	public void close() {
		try {
			queue.put("");
			queue = null;
		} catch (InterruptedException e) {
		}
	}
	
	private class WriterRunnable implements Runnable {
		public void run() {
			boolean shouldWriteHeader = !file.exists();
			try {
				file.getParentFile().mkdirs();
				try (FileWriter writer = new FileWriter(file, true)) {
					if (shouldWriteHeader) {
						writeHeader(writer);
					}
					String str = queue.take();
					while (!str.equals("")) {
						try {
							writer.write(str);
						} catch (Exception e) {
							// consume all exceptions
						}
						str = queue.take();
					}
				}
			} catch (Exception e) {
				log.debug(e);
			}
		}
	}
}
