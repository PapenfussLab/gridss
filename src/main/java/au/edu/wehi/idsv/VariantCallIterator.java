package au.edu.wehi.idsv;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingDeque;
import java.util.concurrent.LinkedBlockingDeque;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Supplier;

import au.edu.wehi.idsv.configuration.VisualisationConfiguration;
import au.edu.wehi.idsv.util.DuplicatingIterable;
import au.edu.wehi.idsv.visualisation.PositionalDeBruijnGraphTracker;
import au.edu.wehi.idsv.visualisation.StateTracker;
import au.edu.wehi.idsv.visualisation.TrackedState;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.lang3.tuple.Pair;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
/**
 * Calls breakpoints from the given evidence
 * 
 * @author Daniel Cameron
 */
public class VariantCallIterator implements CloseableIterator<VariantContextDirectedEvidence> {
	private static final Log log = Log.getInstance(VariantCallIterator.class);
	private static final int ITERATOR_BUFFER_SIZE = 256;
	private final VariantContextDirectedEvidence endOfStream;
	private final ProcessingContext processContext;
	private final DuplicatingIterable<DirectedEvidence> iterable;
	private final QueryInterval[] filterInterval;
	private final BlockingDeque<VariantContextDirectedEvidence> outBuffer = new LinkedBlockingDeque<>(ITERATOR_BUFFER_SIZE);
	private VariantContextDirectedEvidence outBufferHeadNextValidRecord = null;
	private final List<AsyncDirectionalIterator> async = new ArrayList<>();
	private int activeIterators;
	private volatile Exception workerThreadException;
	private VariantCallIterator(ProcessingContext processContext, Iterator<DirectedEvidence> evidence, QueryInterval[] interval, int intervalNumber) {
		this.endOfStream = (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(processContext)
				.id("sentinel")
				.chr(processContext.getReference().getSequenceDictionary().getSequence(0).getSequenceName())
				.start(1)
				.stop(1)
				.alleles("N", "N.")
				.make();
		this.processContext = processContext;
		boolean callBreakends = processContext.getVariantCallingParameters().callBreakends;
		this.activeIterators = callBreakends ? 6 : 4;
		this.iterable = new DuplicatingIterable<>(activeIterators, evidence, ITERATOR_BUFFER_SIZE);
		this.filterInterval = interval;
		for (BreakendDirection localDir : BreakendDirection.values()) {
			for (BreakendDirection remoteDir : BreakendDirection.values()) {
				MaximalEvidenceCliqueIterator it = new MaximalEvidenceCliqueIterator(
						processContext,
						this.iterable.iterator(),
						localDir,
						remoteDir,
						new SequentialIdGenerator(String.format("gridss%d%s%s_", Math.max(intervalNumber, 0), localDir.toChar(), remoteDir.toChar())));
				async.add(new AsyncDirectionalIterator(it, localDir, remoteDir));
			}
			if (callBreakends) {
				BreakendMaximalEvidenceCliqueIterator it = new BreakendMaximalEvidenceCliqueIterator(
						processContext,
						this.iterable.iterator(),
						localDir,
						new SequentialIdGenerator(String.format("gridss%d%s_", Math.max(intervalNumber, 0), localDir.toChar())));
				async.add(new AsyncDirectionalIterator(it, localDir, null));
			}
		}
	}
	public VariantCallIterator(ProcessingContext processContext, Iterator<DirectedEvidence> evidence) {
		this(processContext, evidence, null, -1);
	}
	public VariantCallIterator(AggregateEvidenceSource source) {
		this(source.getContext(), source.iterator(), null, -1);
	}
	public VariantCallIterator(AggregateEvidenceSource source, QueryInterval[] interval, int intervalNumber) {
		this(source.getContext(),
				source.iterator(QueryIntervalUtil.padIntervals(source.getContext().getDictionary(), interval, source.getMaxConcordantFragmentSize() + 1)),
				QueryIntervalUtil.padIntervals(source.getContext().getDictionary(), interval, source.getMaxConcordantFragmentSize() + 1),
				intervalNumber);
	}
	public class AsyncDirectionalIterator<T extends VariantContextDirectedEvidence> implements TrackedState, Closeable {
		private Iterator<T> it;
		private StateTracker currentTracker = null;
		private Collection<TrackedState> currentTrackedObjects = null;
		private T lastElement = null;
		private Thread thread;
		private volatile boolean shouldAbortImmediately = false;
		public AsyncDirectionalIterator(Iterator<T> iterator, BreakendDirection dir1, BreakendDirection dir2) {
			this.it = iterator;
			String positionComponent = (filterInterval == null || filterInterval.length == 0) ? "" : String.format("_%s_%d",
					processContext.getDictionary().getSequence(filterInterval[0].referenceIndex).getSequenceName(),
					filterInterval[0].start);
			if (processContext.getConfig().getVisualisation().maxCliqueTelemetry && it instanceof TrackedState) {
				TrackedState ts = (TrackedState)it;
				String filename = String.format("maxclique%s_%s%s.csv", positionComponent, dir1.toChar(), dir2 == null ? "" : dir2.toChar());
				File file = new File(processContext.getConfig().getVisualisation().directory, filename);
				try {
					this.currentTracker = new StateTracker(file);
					this.currentTrackedObjects = Lists.newArrayList(Iterables.concat(ts.trackedObjects(), this.trackedObjects()));
					this.currentTracker.writeHeader(currentTrackedObjects);
				} catch (IOException e) {
					log.debug("Telemetry failure", e);
				}
			}
			this.it = filterInterval == null ? this.it : wrapFilter(filterInterval, this.it);
			this.thread = new Thread(() -> run());
			this.thread.setDaemon(true);
			this.thread.setName("CallVariants " + positionComponent + dir1.toChar() + (dir2 == null ? "" : dir2.toChar()));
			this.thread.start();
		}
		private Iterator<T> wrapFilter(QueryInterval[] filterInterval, Iterator<T> it) {
			return Iterators.filter(it, v -> {
				if (v instanceof DirectedBreakpoint) {
					BreakpointSummary bs = ((DirectedBreakpoint)v).getBreakendSummary();
					return QueryIntervalUtil.overlaps(filterInterval, bs.referenceIndex, bs.start) ||
							QueryIntervalUtil.overlaps(filterInterval, bs.referenceIndex2, bs.start2);
				} else {
					BreakendSummary be = v.getBreakendSummary();
					return QueryIntervalUtil.overlaps(filterInterval, be.referenceIndex, be.start);
				}
			});
		}
		public void run() {
			try {
				while (it.hasNext() && !shouldAbortImmediately) {
					lastElement = it.next();
					outBuffer.putLast(lastElement);
					if (currentTracker != null) {
						try {
							currentTracker.track(currentTrackedObjects);
						} catch (IOException e) {
							log.debug("Telemetry failure", e);
						}
					}
				}
				outBuffer.putLast(endOfStream);
				if (currentTracker != null) {
					try {
						currentTracker.close();
					} catch (IOException e) {
						log.debug("Telemetry failure during close()", e);
					}
				}
			} catch (Exception e) {
				workerThreadException = e;
			}
		}
		@Override
		public String[] trackedNames() {
			return new String[] {
					"chr",
					"pos",
			};
		}

		@Override
		public Object[] trackedState() {
			return new Object[] {
					lastElement == null ? "" : lastElement.getContig(),
					lastElement == null ? 0 : lastElement.getStart(),
			};
		}

		@Override
		public Collection<TrackedState> trackedObjects() {
			return ImmutableList.of(this);
		}

		@Override
		public void close() {
			shouldAbortImmediately = true;
		}
	}

	private void ensureNext() {
		if (outBufferHeadNextValidRecord == null) {
			try {
				while (activeIterators > 0 && workerThreadException == null) {
					VariantContextDirectedEvidence nextElement = outBuffer.takeFirst();
					if (nextElement != endOfStream) {
						outBufferHeadNextValidRecord = nextElement;
						return;
					} else {
						activeIterators--;
					}
				}
				if (workerThreadException != null) {
					throw new RuntimeException(workerThreadException);
				}
			} catch (InterruptedException e) {
				log.error(e);
				throw new RuntimeException(e);
			}
		}
	}

	@Override
	public boolean hasNext() {
		ensureNext();
		return outBufferHeadNextValidRecord != null;
	}

	@Override
	public VariantContextDirectedEvidence next() {
		ensureNext();
		if (outBufferHeadNextValidRecord == null) {
			throw new NoSuchElementException();
		}
		VariantContextDirectedEvidence result = outBufferHeadNextValidRecord;
		outBufferHeadNextValidRecord = null;
		return result;
	}

	@Override
	public void close() {
		for (AsyncDirectionalIterator adi : async) {
			adi.close();
		}
		CloserUtil.close(iterable);
	}
}
 