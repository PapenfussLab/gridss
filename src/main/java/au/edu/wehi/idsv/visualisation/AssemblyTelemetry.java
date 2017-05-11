package au.edu.wehi.idsv.visualisation;

import java.io.Closeable;
import java.io.File;
import java.io.FileWriter;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import au.edu.wehi.idsv.BreakendDirection;
import htsjdk.samtools.util.Log;

public class AssemblyTelemetry implements Closeable {
	private static final Log log = Log.getInstance(AssemblyTelemetry.class);
	private BlockingQueue<String> queue;
	private File file;
	public AssemblyTelemetry(File telemetryFile) {
		this.file = telemetryFile;
		this.queue = new ArrayBlockingQueue<>(4096);
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
		public void loadGraph(int referenceIndex, int start, int end, int nodes, boolean filtered) {
			String str = String.format("%d,%s,load,%d,%d,%d,%d,%b\n", chunk, direction.toChar(), referenceIndex, start, end, nodes, filtered);
			put(str);
		}

		public void flushContigs(int referenceIndex, int flushStart, int flushEnd, int contigsFlushed) {
			String str = String.format("%d,%s,flushContigs,%d,%d,%d,%d,\n", chunk, direction.toChar(), referenceIndex, flushStart, flushEnd, contigsFlushed);
			put(str);
		}

		public void flushReferenceNodes(int referenceIndex, int flushStart, int flushEnd, int readsFlushed) {
			String str = String.format("%d,%s,flushReferenceNodes,%d,%d,%d,%d,\n", chunk, direction.toChar(), referenceIndex, flushStart, flushEnd, readsFlushed);
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
