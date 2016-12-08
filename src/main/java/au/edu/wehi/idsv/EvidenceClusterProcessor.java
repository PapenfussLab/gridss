package au.edu.wehi.idsv;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Optional;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Executor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
/**
 * Calls breakpoints from the given evidence
 * 
 * @author Daniel Cameron
 */
public class EvidenceClusterProcessor implements Iterator<VariantContextDirectedEvidence> {
	private static final Log log = Log.getInstance(EvidenceClusterProcessor.class);
	private final BlockingQueue<Optional<VariantContextDirectedEvidence>> callBuffer = new ArrayBlockingQueue<Optional<VariantContextDirectedEvidence>>(4096);
	private final ProcessingContext processContext;
	private final AtomicInteger outstandingTasks = new AtomicInteger(0);
	private VariantContextDirectedEvidence nextRecord = null;
	private volatile Exception backgroundThreadException = null;
	public EvidenceClusterProcessor(Executor threadpool, ProcessingContext processContext, Iterable<DirectedEvidence> evidence) throws InterruptedException {
		this.processContext = processContext;
		List<Runnable> tasks = new ArrayList<>(4);
		addToTasksList(tasks, evidence);
		start(threadpool, tasks);
	}
	public EvidenceClusterProcessor(Executor threadpool, AggregateEvidenceSource source) {
		this.processContext = source.getContext();
		List<QueryInterval> chunks = source.getContext().getReference().getIntervals(source.getContext().getConfig().chunkSize);
		List<Runnable> tasks = new ArrayList<>(4 * chunks.size());
		for (int i = 0; i < chunks.size(); i++) {
			addToTasksList(tasks, source, chunks.get(i), i);
		}
		start(threadpool, tasks);
	}
	public void addToTasksList(List<Runnable> list, AggregateEvidenceSource source, QueryInterval interval, int intervalNumber) {
		list.add(() -> { runTask(source, BreakendDirection.Forward, BreakendDirection.Forward, new SequentialIdGenerator(String.format("gridss%dff", intervalNumber)), interval);});
		list.add(() -> { runTask(source, BreakendDirection.Forward, BreakendDirection.Backward, new SequentialIdGenerator(String.format("gridss%dfb", intervalNumber)), interval);});
		list.add(() -> { runTask(source, BreakendDirection.Backward, BreakendDirection.Forward, new SequentialIdGenerator(String.format("gridss%dbf", intervalNumber)), interval);});
		list.add(() -> { runTask(source, BreakendDirection.Backward, BreakendDirection.Backward, new SequentialIdGenerator(String.format("gridss%dbb", intervalNumber)), interval);});
	}
	public void addToTasksList(List<Runnable> list, Iterable<DirectedEvidence> it) {
		list.add(() -> { runTask(it, BreakendDirection.Forward, BreakendDirection.Forward, new SequentialIdGenerator("gridssff"));});
		list.add(() -> { runTask(it, BreakendDirection.Forward, BreakendDirection.Backward, new SequentialIdGenerator("gridssfb"));});
		list.add(() -> { runTask(it, BreakendDirection.Backward, BreakendDirection.Forward, new SequentialIdGenerator("gridssbf"));});
		list.add(() -> { runTask(it, BreakendDirection.Backward, BreakendDirection.Backward, new SequentialIdGenerator("gridssbb"));});
	}
	public void runTask(AggregateEvidenceSource source, BreakendDirection lowDir, BreakendDirection highDir, VariantIdGenerator generator, QueryInterval filter) {
		log.debug("Start processing " + filter);
		try {
			int expandBy = source.getMaxConcordantFragmentSize() + 1;
			try (CloseableIterator<DirectedEvidence> it = source.iterator(new QueryInterval(filter.referenceIndex,  filter.start - expandBy, filter.end + expandBy))) {
				callMaximalCliques(it, lowDir, highDir, generator, filter);
			}
		} catch (Exception e) {
			log.debug(e);
			backgroundThreadException = e;
			throw e;
		} finally {
			outstandingTasks.decrementAndGet();
			callBuffer.add(Optional.empty());
			log.debug("Completed processing " + filter);
		}
	}
	public void runTask(Iterable<DirectedEvidence> it, BreakendDirection lowDir, BreakendDirection highDir, VariantIdGenerator generator) {
		try {
			callMaximalCliques(it.iterator(), lowDir, highDir, generator, null);
		} catch (Exception e) {
			log.debug(e);
			backgroundThreadException = e;
			throw e;
		} finally {
			outstandingTasks.decrementAndGet();
			callBuffer.add(Optional.empty());
		}
	}
	public void callMaximalCliques(Iterator<DirectedEvidence> it, BreakendDirection lowDir, BreakendDirection highDir, VariantIdGenerator generator, QueryInterval filter) {
		MaximalEvidenceCliqueIterator cliqueIt = new MaximalEvidenceCliqueIterator(processContext, it, lowDir, highDir, generator);
		if (filter != null) {
			while (cliqueIt.hasNext()) {
				VariantContextDirectedEvidence v = cliqueIt.next();
				BreakendSummary bs = v.getBreakendSummary();
				if (bs.start >= filter.start && bs.start <= filter.end) {
					callBuffer.add(Optional.of(v));
				}
			}
		} else {
			while (cliqueIt.hasNext()) {
				VariantContextDirectedEvidence v = cliqueIt.next();
				callBuffer.add(Optional.of(v));
			}
		}
	}
	private void start(Executor threadpool, List<Runnable> tasks) {
		if (threadpool == null) throw new NullPointerException("threadpool cannot be null");
		assert(nextRecord == null);
		assert(callBuffer.isEmpty());
		outstandingTasks.set(tasks.size());
		for (Runnable task : tasks) {
			threadpool.execute(task);
		}
	}
	private void ensureNextRecord() throws InterruptedException {
		
		if (nextRecord != null) return;
		// If there's still a task running then we know at least one
		// worker thread will write a null sentinal
		while (outstandingTasks.get() > 0) {
			Optional<VariantContextDirectedEvidence> result = callBuffer.poll(1, TimeUnit.SECONDS);
			if (backgroundThreadException != null) return;
			if (result != null && result.isPresent()) {
				nextRecord = result.get();
				return;
			}
		}
		// If the workers are all done, then we're never going to
		// get any more records added to the queue
		// (actually, technically we can, but they're sentinal records which we can ignore)
		while (!callBuffer.isEmpty()) {
			Optional<VariantContextDirectedEvidence> result = callBuffer.poll();
			if (result.isPresent()) {
				nextRecord = result.get();
				return;
			}
		}
	}
	@Override
	public boolean hasNext() {
		try {
			ensureNextRecord();
			if (backgroundThreadException != null) {
				throw new RuntimeException(backgroundThreadException);
			}
		} catch (InterruptedException e) {
			log.error(e, "Unexpectedly interrupted waiting for maximal cliques to be calculated.");
			throw new RuntimeException(e);
		}
		return nextRecord != null;
	}
	@Override
	public VariantContextDirectedEvidence next() {
		if (!hasNext()) throw new NoSuchElementException();
		VariantContextDirectedEvidence result = nextRecord;
		nextRecord = null;
		return result;
	}
}
 