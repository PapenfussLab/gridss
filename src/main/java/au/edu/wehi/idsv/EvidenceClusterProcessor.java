package au.edu.wehi.idsv;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;

import java.util.Iterator;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;

import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.util.DuplicatingIterable;

import com.google.common.collect.AbstractIterator;
/**
 * Calls breakpoints from the given evidence
 * 
 * @author Daniel Cameron
 */
public class EvidenceClusterProcessor extends AbstractIterator<VariantContextDirectedEvidence> implements CloseableIterator<VariantContextDirectedEvidence> {
	private static final Log log = Log.getInstance(EvidenceClusterProcessor.class);
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
		DuplicatingIterable<DirectedEvidence> dib = new DuplicatingIterable<DirectedEvidence>(4, this.underlying, 32);
		this.threads = new MaximalCliqueIteratorRunnable[] {
			new MaximalCliqueIteratorRunnable(context, dib.iterator(), BreakendDirection.Forward, BreakendDirection.Forward),
			new MaximalCliqueIteratorRunnable(context, dib.iterator(), BreakendDirection.Forward, BreakendDirection.Backward),
			new MaximalCliqueIteratorRunnable(context, dib.iterator(), BreakendDirection.Backward, BreakendDirection.Forward),
			new MaximalCliqueIteratorRunnable(context, dib.iterator(), BreakendDirection.Backward, BreakendDirection.Backward)
		};
		for (MaximalCliqueIteratorRunnable r : threads) {
			r.setUncaughtExceptionHandler(handler);
			r.start();
		}
	}
	private class MaximalCliqueIteratorRunnable extends Thread {
		private final MaximalEvidenceCliqueIterator cliqueIt;
		public MaximalCliqueIteratorRunnable(ProcessingContext processContext, Iterator<DirectedEvidence> evidenceIt, BreakendDirection lowDir, BreakendDirection highDir) {
			cliqueIt = new MaximalEvidenceCliqueIterator(processContext, evidenceIt, lowDir, highDir);
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
 