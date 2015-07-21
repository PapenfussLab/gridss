package au.edu.wehi.idsv.visualisation;

import htsjdk.samtools.util.CloserUtil;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.lang.ref.WeakReference;
import java.nio.charset.StandardCharsets;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.List;

import au.edu.wehi.idsv.visualisation.TrackedBuffer.NamedTrackedBuffer;

public class BufferTracker {
	private final List<WeakReference<TrackedBuffer>> bufferObjects = new ArrayList<WeakReference<TrackedBuffer>>();
	private final File output;
	private final int writeIntervalInSeconds;
	public BufferTracker(File output, int writeIntervalInSeconds) {
		this.output = output;
		this.writeIntervalInSeconds = writeIntervalInSeconds;
		Worker worker = new Worker();
		worker.setName("BufferTracker");
		worker.setDaemon(false);
		worker.start();
	}
	public synchronized void register(String context, TrackedBuffer obj) {
		obj.setTrackedBufferContext(context);
		bufferObjects.add(new WeakReference<TrackedBuffer>(obj));
	}
	private synchronized String getCsvRows() {
		StringBuilder sb = new StringBuilder();
		String timestamp = LocalDateTime.now().toString();
		for (WeakReference<TrackedBuffer> wr : bufferObjects) {
			TrackedBuffer buffer = wr.get();
			if (buffer != null) {
				for (NamedTrackedBuffer bufferSize : buffer.currentTrackedBufferSizes()) {
					sb.append(timestamp);
					sb.append(',');
					sb.append(bufferSize.name);
					sb.append(',');
					sb.append(Long.toString(bufferSize.size));
					sb.append('\n');
				}
			}
		}
		return sb.toString();
	}
	private void append() {
		FileOutputStream os = null;
		String str = getCsvRows();
		if (!str.isEmpty()) {
			try {
				os = new FileOutputStream(output, true);
				os.write(str.getBytes(StandardCharsets.UTF_8));
				os.flush();
				os.close();
				os = null;
			} catch (IOException e) {
			} finally {
				CloserUtil.close(os);
			}
		}
	}
	private class Worker extends Thread {
		@Override
		public void run() {
			while (true) {
				try {
					Thread.sleep(writeIntervalInSeconds * 1000);
					append();
				} catch (InterruptedException e) {
					// swallow
				}
			}
		}
	}
}
