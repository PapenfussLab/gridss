package au.edu.wehi.idsv;

import java.util.Iterator;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import com.google.common.collect.AbstractIterator;

import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.util.DuplicatingIterable;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
/**
 * Calls breakpoints from the given evidence
 * 
 * @author Daniel Cameron
 */
public class EvidenceClusterProcessor extends AbstractIterator<VariantContextDirectedEvidence> implements CloseableIterator<VariantContextDirectedEvidence> {
	private static final Log log = Log.getInstance(EvidenceClusterProcessor.class);
	private static final int EVIDENCE_BUFFER_SIZE = 1024;
	private static final AtomicInteger threadCount = new AtomicInteger(0);
	private final BlockingQueue<VariantContextDirectedEvidence> callBuffer = new ArrayBlockingQueue<VariantContextDirectedEvidence>(64);
	private MaximalCliqueIteratorRunnable[] threads;
	private AutoClosingIterator<DirectedEvidence> underlying;
	private volatile Throwable uncaught = null;
	private volatile boolean isClosed = false;
	private Thread.UncaughtExceptionHandler handler = new Thread.UncaughtExceptionHandler() {
	    public void uncaughtException(Thread th, Throwable ex) {
	    	log.error(ex);
	    	uncaught = ex;
	    }
	};
	public EvidenceClusterProcessor(ProcessingContext context, Iterator<DirectedEvidence> evidence) {
		if (evidence == null) throw new IllegalArgumentException();
		this.underlying = new AutoClosingIterator<DirectedEvidence>(evidence);
		// Start each on their own thread
		DuplicatingIterable<DirectedEvidence> dib = new DuplicatingIterable<DirectedEvidence>(4, this.underlying, EVIDENCE_BUFFER_SIZE);
		this.threads = new MaximalCliqueIteratorRunnable[] {
			new MaximalCliqueIteratorRunnable(context, dib.iterator(), BreakendDirection.Forward, BreakendDirection.Forward, new SequentialIdGenerator("gridssff")),
			new MaximalCliqueIteratorRunnable(context, dib.iterator(), BreakendDirection.Forward, BreakendDirection.Backward, new SequentialIdGenerator("gridssfb")),
			new MaximalCliqueIteratorRunnable(context, dib.iterator(), BreakendDirection.Backward, BreakendDirection.Forward, new SequentialIdGenerator("gridssbf")),
			new MaximalCliqueIteratorRunnable(context, dib.iterator(), BreakendDirection.Backward, BreakendDirection.Backward, new SequentialIdGenerator("gridssbb"))
		};
		for (MaximalCliqueIteratorRunnable r : threads) {
			r.setUncaughtExceptionHandler(handler);
			r.start();
		}
	}
	
	private class MaximalCliqueIteratorRunnable extends Thread {
		private final MaximalEvidenceCliqueIterator cliqueIt;
		public MaximalCliqueIteratorRunnable(ProcessingContext processContext, Iterator<DirectedEvidence> evidenceIt, BreakendDirection lowDir, BreakendDirection highDir, VariantIdGenerator generator) {
			cliqueIt = new MaximalEvidenceCliqueIterator(processContext, evidenceIt, lowDir, highDir, generator);
			this.setName(String.format("MaxClique%s%s-%d", lowDir.toChar(), highDir.toChar(), threadCount.incrementAndGet()));
		}
		@Override
		public void run() {
			while (cliqueIt.hasNext() && !isClosed) {
				try {
					callBuffer.put(cliqueIt.next());
				} catch (InterruptedException e) {
					if (!isClosed) {
						log.error("Interrupted waiting on output buffer");
						throw new RuntimeException(e);
					}
				}
			}
		}
	}
	private boolean workersDone() {
		boolean allClosed = true;
		for (MaximalCliqueIteratorRunnable r : threads) {
			allClosed &= !r.isAlive();
		}
		return allClosed;
	}
	private void checkBackgroundThreadSuccess() {
		if (uncaught != null) {
			close();
			throw new RuntimeException(uncaught);
		}
	}
	@Override
	protected VariantContextDirectedEvidence computeNext() {
		checkBackgroundThreadSuccess();
		if (isClosed) return endOfData();
		// TODO: improve on busy wait loop
		while (!workersDone()) {
			VariantContextDirectedEvidence result;
			try {
				result = callBuffer.poll(500, TimeUnit.MILLISECONDS);
				checkBackgroundThreadSuccess();
				if (result != null) {
					return result;
				}
			} catch (InterruptedException e) {
				log.debug("Interrupted waiting for maximal clique call");
				checkBackgroundThreadSuccess();
			}
		}
		// flush remaining buffer
		while (!callBuffer.isEmpty()) return callBuffer.poll();
		close();
		return endOfData();
	}
	@Override
	public void close() {	
		isClosed = true;
		if (underlying != null) {
			underlying.close();
			underlying = null;
		}
		callBuffer.clear();
		if (threads != null) {
			for (int i = 0; i < threads.length; i++) {
				try {
					if (threads[i] != null) {
						threads[i].join(1000);
					}
				} catch (InterruptedException e) {
					log.warn("Interrupted waiting for child thread to join");
				}
				threads[i] = null;
			}
			threads = null;
		}
	}
}
 