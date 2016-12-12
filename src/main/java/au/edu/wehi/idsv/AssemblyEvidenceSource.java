package au.edu.wehi.idsv;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Optional;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import com.google.common.base.Stopwatch;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.util.concurrent.MoreExecutors;

import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.configuration.AssemblyConfiguration;
import au.edu.wehi.idsv.debruijn.positional.PositionalAssembler;
import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.sam.SAMFileUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.FileHelper;
import gridss.SoftClipsToSplitReads;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

/**
 * Structural variant supporting contigs generated from assembly
 * 
 *  TODO: define separate functions for the attributes of the assemblies, and the source reads.
 *  It is currently unclear whether getMaxReadLength() refers to the assembly or source
 *  reads (it is the latter).
 *  
 * @author Daniel Cameron
 *
 */
public class AssemblyEvidenceSource extends SAMEvidenceSource {
	private static final Log log = Log.getInstance(AssemblyEvidenceSource.class);
	private final List<SAMEvidenceSource> source;
	private final int maxSourceFragSize;
	private final int minSourceFragSize;
	private final int maxReadLength;
	private final int maxMappedReadLength;
	private final IntervalBed throttled;
	/**
	 * Generates assembly evidence based on the given evidence
	 * @param evidence evidence for creating assembly
	 * @param intermediateFileLocation location to store intermediate files
	 */
	public AssemblyEvidenceSource(ProcessingContext processContext, List<SAMEvidenceSource> evidence, File assemblyFile) {
		super(processContext, assemblyFile, null, -1);
		this.source = evidence;
		this.maxSourceFragSize = evidence.stream().mapToInt(s -> s.getMaxConcordantFragmentSize()).max().orElse(0);
		this.minSourceFragSize = evidence.stream().mapToInt(s -> s.getMinConcordantFragmentSize()).min().orElse(0);
		this.maxReadLength = evidence.stream().mapToInt(s -> s.getMaxReadLength()).max().orElse(0);
		this.maxMappedReadLength = evidence.stream().mapToInt(s -> Math.max(s.getMaxReadLength(), s.getMaxReadMappedLength())).max().orElse(0);
		assert(maxSourceFragSize >= minSourceFragSize);
		this.throttled = new IntervalBed(getContext().getDictionary(), getContext().getLinear());
	}
	/**
	 * Perform breakend assembly 
	 * @param threadpool 
	 * @throws IOException 
	 */
	public void assembleBreakends(ExecutorService threadpool) throws IOException {
		if (threadpool == null) {
			threadpool = MoreExecutors.newDirectExecutorService();
		}
		SAMFileHeader header = getContext().getBasicSamHeader();
		header.setSortOrder(SortOrder.unsorted);
		// TODO: add assembly @PG header
		File tmpout = FileSystemContext.getWorkingFileFor(getFile());
		File filteredout = FileSystemContext.getWorkingFileFor(getFile(), "filtered.");
		AssemblyIterator it = new AssemblyIterator(threadpool);
		try (SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, tmpout)) {
			try (SAMFileWriter filteredWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, filteredout)) {
				while (it.hasNext()) {
					SAMRecord asm = it.next();
					asm = transformAssembly(asm);
					if (shouldFilterAssembly(asm)) {
						filteredWriter.addAlignment(asm);
					} else {
						writer.addAlignment(asm);
					}
				}
			}
	    }
		SAMFileUtil.sort(getContext().getFileSystemContext(), tmpout, getFile(), SortOrder.coordinate);
		if (gridss.Defaults.DELETE_TEMPORARY_FILES) {
			FileHelper.delete(tmpout, true);
		}
		if (!getContext().getAssemblyParameters().writeFiltered) {
			FileHelper.delete(filteredout, true);
		}
		if (filteredout.exists()) {
			File filteredsorted = FileSystemContext.getWorkingFileFor(getFile(), "filtered.sorted.");
			SAMFileUtil.sort(getContext().getFileSystemContext(), filteredout, filteredsorted, SortOrder.coordinate);
			FileHelper.move(filteredsorted, filteredout, true);
		}
		File throttledFilename = new File(getFile().getAbsolutePath() + ".throttled.bed");
		try {
			if (throttled.size() > 0) {
				throttled.write(throttledFilename, "Regions of high coverage where only a portion of supporting reads considered for assembly");
			}
		} catch (IOException e) {
			log.warn(e, "Unable to write " + throttledFilename.getAbsolutePath());
		}
	}
	@Override
	public void ensureExtracted() throws IOException {
		ensureMetrics();
		File svFile = getContext().getFileSystemContext().getSVBam(getFile());
		File withsplitreadsFile = FileSystemContext.getWorkingFileFor(svFile, "gridss.tmp.withsplitreads.");
		ensureMetrics();
		// Regenerate from from the intermediate file furtherest through the pipeline
		// extract -> query sort -> tag -> split read -> back to coordinate sorted
		// We want to tag before generating split reads so all splits are guaranteed to
		// have the same tags
		if (!svFile.exists()) {
			log.info("Identifying split reads for " + getFile().getAbsolutePath());
			List<String> args = Lists.newArrayList(
					"INPUT=" + getFile().getAbsolutePath(),
					"OUTPUT=" + svFile.getAbsolutePath());
			execute(new SoftClipsToSplitReads(), args);
		}
		SAMFileUtil.sort(getContext().getFileSystemContext(), withsplitreadsFile, svFile, SortOrder.coordinate);
	}
	public boolean shouldFilterAssembly(SAMRecord asm) {
		AssemblyConfiguration ap = getContext().getAssemblyParameters();
		AssemblyAttributes attr = new AssemblyAttributes(asm);
		// reference assembly
		List<SingleReadEvidence> breakends = SingleReadEvidence.createEvidence(this, asm);
		if (breakends.size() == 0) {
			return true;
		}
		// too long
		int breakendLength = SAMRecordUtil.getSoftClipLength(asm, attr.getAssemblyDirection());
		if (breakendLength > ap.maxExpectedBreakendLengthMultiple * getMaxConcordantFragmentSize()) {
			log.debug(String.format("Filtering %s at %s:%d due to misassembly (breakend %dbp)",
					asm.getReadName(),
					asm.getReferenceName(), asm.getAlignmentStart(),
					breakendLength));
			return true;
		}		
		// too few reads
		if (attr.getAssemblySupportCount() < ap.minReads) {
			return true;
		}
		// unanchored assembly that not actually any longer than any of the reads that were assembled together
		if (attr.getAssemblySupportCountSoftClip() == 0 && breakendLength <= attr.getAssemblyReadPairLengthMax()) {
			// assembly length = 1 read
			// at best, we've just error corrected a single reads with other reads.
			// at worst, we've created a misassembly.
			return true;
		}
		return false;
	}
	public SAMRecord transformAssembly(SAMRecord assembly) {
		if (assembly.getReadUnmappedFlag()) return assembly;
		if (CigarUtil.widthOfImprecision(assembly.getCigar()) == 0) {
			// some assemblies actually match the reference and we can ignore these
			// such assemblies are caused by sequencing errors or SNVs causing
			// spurious soft clips
			SAMRecordUtil.unclipExactReferenceMatches(getContext().getReference(), assembly);
			/*
			float anchorMatchRate = SAMRecordUtil.getAlignedIdentity(assembly);
			if (anchorMatchRate < getContext().getAssemblyParameters().minAnchorIdentity) {
				SAMRecord realigned = SAMRecordUtil.realign(getContext().getReference(), assembly, getContext().getConfig().getAssembly().realignmentWindowSize, true);
				// TODO: check if realignment is actually better
				// (don't allow very short anchors)
				// (don't allow flipping anchor to other side)
				// ...
				float realignedAnchorMatchRate = SAMRecordUtil.getAlignedIdentity(realigned);
				if (realignedAnchorMatchRate >= getContext().getAssemblyParameters().minAnchorIdentity &&
						CigarUtil.commonReferenceBases(realigned.getCigar(), assembly.getCigar()) >= CigarUtil.countMappedBases(assembly.getCigar().getCigarElements()) * getContext().getConfig().getAssembly().minPortionOfAnchorRetained) {
					assembly = realigned;
				}
			}
			*/
		}
		return assembly;
	}
	private Iterator<DirectedEvidence> throttled(Iterator<DirectedEvidence> it) {
		AssemblyConfiguration ap = getContext().getAssemblyParameters();
		DirectedEvidenceDensityThrottlingIterator dit = new DirectedEvidenceDensityThrottlingIterator(
				throttled,
				getContext().getDictionary(),
				getContext().getLinear(),
				it,
				Math.max(ap.downsampling.minimumDensityWindowSize, getMaxConcordantFragmentSize()),
				ap.downsampling.acceptDensityPortion * ap.downsampling.targetEvidenceDensity,
				ap.downsampling.targetEvidenceDensity);
		getContext().registerBuffer(AssemblyEvidenceSource.class.getName() + ".throttle", dit);
		return dit;
	}
	public int getMaxAssemblyLength() {
		// We could extract by generating metrics for the assembly
		float maxExpected = getContext().getAssemblyParameters().maxExpectedBreakendLengthMultiple * getMaxConcordantFragmentSize();
		maxExpected = Math.max(maxExpected, getContext().getAssemblyParameters().anchorLength);
		maxExpected *= 2; // anchor + breakend length
		return (int)maxExpected;
	}
	@Override
	public int getMaxConcordantFragmentSize() {
		return maxSourceFragSize;
	}
	@Override
	public int getMinConcordantFragmentSize() {
		return minSourceFragSize;
	}
	/**
	 * Maximum read length of reads contributing to assemblies
	 * @return
	 */
	public int getMaxReadLength() {
		return maxReadLength;
	}
	public int getMaxMappedReadLength() {
		return maxMappedReadLength;
	}
	@SuppressWarnings("unused")
	private Iterator<SAMRecord> getAllAssemblies_single_threaded() {
		ProgressLogger progressLog = new ProgressLogger(log);
		List<Iterator<SAMRecord>> list = new ArrayList<>();
		for (BreakendDirection direction : BreakendDirection.values()) {
			CloseableIterator<DirectedEvidence> it = mergedIterator(source, false);
			Iterator<DirectedEvidence> throttledIt = throttled(it);
			ProgressLoggingDirectedEvidenceIterator<DirectedEvidence> loggedIt = new ProgressLoggingDirectedEvidenceIterator<>(getContext(), throttledIt, progressLog);
			Iterator<SAMRecord> evidenceIt = new PositionalAssembler(getContext(), this, loggedIt, direction);
	    	list.add(evidenceIt);
		}
		return Iterators.concat(list.iterator());
	}
	private class AssemblyIterator implements Iterator<SAMRecord> {
		private final LinkedBlockingQueue<Optional<SAMRecord>> outputQueue = new LinkedBlockingQueue<>();
		private final AtomicInteger outstandingTasks;
		private SAMRecord nextRecord = null;
		private volatile Exception backgroundThreadException = null;
		public AssemblyIterator(ExecutorService threadpool) {
			if (threadpool == null) throw new NullPointerException("threadpool cannot be null");
			List<QueryInterval> chunks = getContext().getReference().getIntervals(getContext().getConfig().chunkSize);
			this.outstandingTasks = new AtomicInteger(2 * chunks.size()); // a forward and backward for each chunk 
			assembleTasks(threadpool, chunks);
		}
		private void assembleTasks(ExecutorService threadpool, List<QueryInterval> chunks) {
			for (QueryInterval currentChunck : chunks) {
				final QueryInterval qi = currentChunck;
				threadpool.submit(() -> { assembleChunk(qi, BreakendDirection.Forward); });
				threadpool.submit(() -> { assembleChunk(qi, BreakendDirection.Backward); });
			}
		}
		private void assembleChunk(QueryInterval qi, BreakendDirection direction) {
			String chuckName = String.format("%s:%d-%d %s", getContext().getReference().getSequenceDictionary().getSequence(qi.referenceIndex).getSequenceName(), qi.start, qi.end, direction);
			try {
				log.info(String.format("Starting assembly on interval %s", chuckName));
				Stopwatch timer = Stopwatch.createStarted();
				int expansion = (int)(2 * getMaxConcordantFragmentSize() * getContext().getConfig().getAssembly().maxExpectedBreakendLengthMultiple) + 1;
				int end = getContext().getReference().getSequenceDictionary().getSequence(qi.referenceIndex).getSequenceLength();
				QueryInterval expanded = new QueryInterval(qi.referenceIndex, Math.max(1, qi.start - expansion), Math.min(end, qi.end + expansion));				
				try (CloseableIterator<DirectedEvidence> input = mergedIterator(source, expanded)) {
					Iterator<DirectedEvidence> throttledIt = throttled(input);
					Iterator<SAMRecord> evidenceIt = new PositionalAssembler(getContext(), AssemblyEvidenceSource.this, throttledIt, direction);
					while (evidenceIt.hasNext()) {
						SAMRecord ass = evidenceIt.next();
						if (ass.getAlignmentStart() >= qi.start && ass.getAlignmentStart() <= qi.end) {
							// only output assemblies that start within our chuck
							outputQueue.put(Optional.of(ass));
						}
					}
					outstandingTasks.decrementAndGet();
					// Sentinal value to indicate we are done
					outputQueue.put(Optional.empty());
				}
				timer.stop();
				log.info(String.format("Completed assembly on interval %s in %ds (%s)", chuckName, timer.elapsed(TimeUnit.SECONDS), timer.toString()));
			} catch (Exception e) {
				log.error(e, "Error assembling ", chuckName);
				backgroundThreadException = e;
				if (getContext().getConfig().terminateOnFirstError) {
					System.exit(1);
				}
			}
		}
		private void ensureNextRecord() throws InterruptedException {
			
			if (nextRecord != null) return;
			// If there's still a task running then we know at least one
			// worker thread will write a null sentinal
			while (outstandingTasks.get() > 0) {
				Optional<SAMRecord> result = outputQueue.poll(1, TimeUnit.SECONDS);
				if (backgroundThreadException != null) return;
				if (result != null && result.isPresent()) {
					nextRecord = result.get();
					return;
				}
			}
			// If the workers are all done, then we're never going to
			// get any more records added to the queue
			// (actually, technically we can, but they're sentinal records which we can ignore)
			while (!outputQueue.isEmpty()) {
				Optional<SAMRecord> result = outputQueue.poll();
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
				log.error(e, "Unexpectedly interrupted waiting for assembly contigs to be generated.");
				throw new RuntimeException(e);
			}
			return nextRecord != null;
		}
		@Override
		public SAMRecord next() {
			if (!hasNext()) throw new NoSuchElementException();
			SAMRecord result = nextRecord;
			nextRecord = null;
			return result;
		}
		
	}
}
