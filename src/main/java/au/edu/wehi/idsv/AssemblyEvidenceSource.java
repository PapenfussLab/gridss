package au.edu.wehi.idsv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.configuration.AssemblyConfiguration;
import au.edu.wehi.idsv.debruijn.positional.DirectedPositionalAssembler;
import au.edu.wehi.idsv.debruijn.positional.PositionalAssembler;
import au.edu.wehi.idsv.debruijn.subgraph.DeBruijnSubgraphAssembler;
import au.edu.wehi.idsv.pipeline.CreateAssemblyReadPair;
import au.edu.wehi.idsv.sam.SAMFileUtil.SortCallable;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.util.AutoClosingMergedIterator;
import au.edu.wehi.idsv.util.FileHelper;
import au.edu.wehi.idsv.validation.OrderAssertingIterator;
import au.edu.wehi.idsv.validation.PairedEvidenceTracker;

import com.google.common.base.Function;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;


public class AssemblyEvidenceSource extends EvidenceSource {
	private static final Log log = Log.getInstance(AssemblyEvidenceSource.class);
	private static final int PROGRESS_UPDATE_INTERVAL_SECONDS = 60;
	private final List<SAMEvidenceSource> source;
	private final int maxSourceFragSize;
	private final int minSourceFragSize;
	private final int maxReadLength;
	private final int maxMappedReadLength;
	private final FileSystemContext fsc;
	/**
	 * Generates assembly evidence based on the given evidence
	 * @param evidence evidence for creating assembly
	 * @param intermediateFileLocation location to store intermediate files
	 */
	public AssemblyEvidenceSource(ProcessingContext processContext, List<SAMEvidenceSource> evidence, File intermediateFileLocation) {
		super(processContext, intermediateFileLocation);
		this.fsc = getContext().getFileSystemContext();
		this.source = evidence;
		this.maxSourceFragSize = evidence.stream().mapToInt(s -> s.getMaxConcordantFragmentSize()).max().orElse(0);
		this.minSourceFragSize = evidence.stream().mapToInt(s -> s.getMinConcordantFragmentSize()).min().orElse(0);
		this.maxReadLength = evidence.stream().mapToInt(s -> s.getMaxReadLength()).max().orElse(0);
		this.maxMappedReadLength = evidence.stream().mapToInt(s -> Math.max(s.getMaxReadLength(), s.getMaxReadMappedLength())).max().orElse(0);
		assert(maxSourceFragSize >= minSourceFragSize);
	}
	@Override
	public int getRealignmentIterationCount() {
		return getContext().getRealignmentParameters().assemblyIterations;
	}
	public void ensureAssembled() {
		ensureAssembled(null);
	}
	/**
	 * Assembles evidence
	 * @param threadpool thread pool containing desired number of active threads
	 */
	public void ensureAssembled(ExecutorService threadpool) {
		if (!isProcessingComplete()) {
			process(threadpool);
		}
		if (isRealignmentComplete()) {
			CreateAssemblyReadPair step = new CreateAssemblyReadPair(getContext(), this, source);
			step.process(EnumSet.of(ProcessStep.SORT_REALIGNED_ASSEMBLIES), threadpool);
			step.close();
		}
	}
	public CloseableIterator<SAMRecordAssemblyEvidence> iterator(final boolean includeRemote, final boolean includeFiltered) {
		CloseableIterator<SAMRecordAssemblyEvidence> it;
		if (getContext().shouldProcessPerChromosome()) {
			// Lazily iterator over each input
			it = new PerChromosomeAggregateIterator<SAMRecordAssemblyEvidence>(getContext().getReference().getSequenceDictionary(), new Function<String, Iterator<SAMRecordAssemblyEvidence>>() {
				@Override
				public Iterator<SAMRecordAssemblyEvidence> apply(String chr) {
					return perChrIterator(includeRemote, includeFiltered, chr);
				}
			});
		} else {
			it = singleFileIterator(includeRemote, includeFiltered);
		}
		if (Defaults.SANITY_CHECK_ITERATORS && includeRemote) {
			it = new PairedEvidenceTracker<SAMRecordAssemblyEvidence>(it);
		}
		return it;
	}
	public CloseableIterator<SAMRecordAssemblyEvidence> iterator(boolean includeRemote, boolean includeFiltered, String chr) {
		if (getContext().shouldProcessPerChromosome()) {
			return perChrIterator(includeRemote, includeFiltered, chr);
		} else {
			return new ChromosomeFilteringIterator<SAMRecordAssemblyEvidence>(singleFileIterator(includeRemote, includeFiltered), getContext().getDictionary().getSequence(chr).getSequenceIndex(), true);
		}
	}
	protected CloseableIterator<SAMRecordAssemblyEvidence> perChrIterator(boolean includeRemote, boolean includeFiltered, String chr) {
		FileSystemContext fsc = getContext().getFileSystemContext();
		if (!isReadPairingComplete()) {
			return samAssemblyRealignIterator(
				includeRemote,
				includeFiltered,
				fsc.getAssemblyRawBamForChr(input, chr),
				IntStream.range(0, getRealignmentIterationCount()).mapToObj(i -> fsc.getRealignmentBamForChr(input, chr, i)).collect(Collectors.toList()),
				chr);
		} else {
			return samReadPairIterator(
					includeRemote,
					includeFiltered,
					fsc.getAssemblyForChr(input, chr),
					fsc.getAssemblyMateForChr(input, chr),
					chr);
		}
	}
	protected CloseableIterator<SAMRecordAssemblyEvidence> singleFileIterator(boolean includeRemote, boolean includeFiltered) {
		FileSystemContext fsc = getContext().getFileSystemContext();
		if (!isReadPairingComplete()) {
			return samAssemblyRealignIterator(
					includeRemote,
					includeFiltered,
					fsc.getAssemblyRawBam(input),
					IntStream.range(0, getRealignmentIterationCount()).mapToObj(i -> fsc.getRealignmentBam(input, i)).collect(Collectors.toList()),
					"");
		} else {
			return samReadPairIterator(
					includeRemote,
					includeFiltered,
					fsc.getAssembly(input),
					fsc.getAssemblyMate(input),
					"");
		}
	}
	private CloseableIterator<SAMRecordAssemblyEvidence> samReadPairIterator(
			boolean includeRemote,
			boolean includeFiltered,
			File assembly,
			File mate,
			String chr) {
		CloseableIterator<SAMRecord> it = getContext().getSamReaderIterator(assembly);
		CloseableIterator<SAMRecord> mateIt = getContext().getSamReaderIterator(mate);
		CloseableIterator<SAMRecordAssemblyEvidence> evidenceIt = new SAMRecordAssemblyEvidenceReadPairIterator(getContext(), this, it, mateIt, includeRemote);
		getContext().registerBuffer("assembly.rp." + chr, (SAMRecordAssemblyEvidenceReadPairIterator)evidenceIt);
		Iterator<SAMRecordAssemblyEvidence> filteredIt = includeFiltered ? evidenceIt : new SAMRecordAssemblyEvidenceFilteringIterator(getContext(), evidenceIt);
		DirectEvidenceWindowedSortingIterator<SAMRecordAssemblyEvidence> dit = new DirectEvidenceWindowedSortingIterator<SAMRecordAssemblyEvidence>(
				getContext(),
				getAssemblyWindowSize(),
				filteredIt);
		getContext().registerBuffer("assembly.rp." + chr, dit);
		CloseableIterator<SAMRecordAssemblyEvidence> sortedIt = new AutoClosingIterator<SAMRecordAssemblyEvidence>(dit, ImmutableList.<Closeable>of(it, mateIt, evidenceIt));
		if (Defaults.SANITY_CHECK_ITERATORS) {
			sortedIt = new AutoClosingIterator<SAMRecordAssemblyEvidence>(new OrderAssertingIterator<SAMRecordAssemblyEvidence>(sortedIt, DirectedEvidenceOrder.ByNatural), ImmutableList.<Closeable>of(sortedIt));
		}
		return sortedIt;
	}
	private CloseableIterator<SAMRecordAssemblyEvidence> samAssemblyRealignIterator(
			boolean includeRemote,
			boolean includeFiltered,
			File breakend,
			List<File> realignment,
			String chr) {
		if (!isProcessingComplete()) {
			log.error("Assemblies not yet generated.");
			throw new IllegalStateException("Assemblies not yet generated");
		}
		if (includeRemote) {
			log.error("Realignment sorting not complete.");
			throw new IllegalStateException("Remote assembly iteration requires realignment complete");
		}
		List<Closeable> toClose = Lists.newArrayList();
		List<CloseableIterator<SAMRecord>> realignedIt = new ArrayList<CloseableIterator<SAMRecord>>();
		if (isRealignmentComplete()) {
			for (File realignmenti : realignment) {
				CloseableIterator<SAMRecord> realignmentiit = getContext().getSamReaderIterator(realignmenti); 
				realignedIt.add(realignmentiit);
				toClose.add(realignmentiit);
			}
		} else {
			log.debug(String.format("Assembly realignment for %s not completed", breakend));
		}
		CloseableIterator<SAMRecord> rawReaderIt = getContext().getSamReaderIterator(breakend, SortOrder.coordinate);
		toClose.add(rawReaderIt);
		CloseableIterator<SAMRecordAssemblyEvidence> evidenceIt = new SAMRecordAssemblyEvidenceIterator(
				getContext(), this,
				new AutoClosingIterator<SAMRecord>(rawReaderIt, realignedIt),
				realignedIt,
				true);
		toClose.add(evidenceIt);
		Iterator<SAMRecordAssemblyEvidence> filteredIt = includeFiltered ? evidenceIt : new SAMRecordAssemblyEvidenceFilteringIterator(getContext(), evidenceIt);
		// Change sort order to breakend position order
		DirectEvidenceWindowedSortingIterator<SAMRecordAssemblyEvidence> dit = new DirectEvidenceWindowedSortingIterator<SAMRecordAssemblyEvidence>(
				getContext(),
				// double window size due to spanning assembly edge case: 
				//   |--- window ---|
				// <-A              A->
				//                  B->              <-B
				//                  ^
				// SAMRecord anchored at breakend
				// since our SAMRecordAssemblyEvidenceIterator returns 4 records (both breakends, both assemblies)
				// the breakend position can be +- window size size away from the start position of the underlying SAMRecord
				2 * getAssemblyWindowSize(),
				filteredIt);
		getContext().registerBuffer("assembly.sam." + chr, dit);
		CloseableIterator<SAMRecordAssemblyEvidence> sortedIt = new AutoClosingIterator<SAMRecordAssemblyEvidence>(dit, toClose);
		if (Defaults.SANITY_CHECK_ITERATORS) {
			sortedIt = new AutoClosingIterator<SAMRecordAssemblyEvidence>(new OrderAssertingIterator<SAMRecordAssemblyEvidence>(sortedIt, DirectedEvidenceOrder.ByNatural), ImmutableList.<Closeable>of(sortedIt));
		}
		return sortedIt;
	}
	private boolean isProcessingComplete() {
		List<File> target = new ArrayList<File>();
		List<File> source = new ArrayList<File>();
		if (getContext().shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : getContext().getReference().getSequenceDictionary().getSequences()) {
				source.add(null); target.add(fsc.getAssemblyRawBamForChr(input, seq.getSequenceName()));
				source.add(null); target.add(fsc.getRealignmentFastqForChr(input, seq.getSequenceName(), 0));
			}
		} else {
			source.add(null); target.add(fsc.getAssemblyRawBam(input));
			source.add(null); target.add(fsc.getRealignmentFastq(input, 0));
		}
		return IntermediateFileUtil.checkIntermediate(target, source, getContext().getConfig().ignoreFileTimestamps);
	}
	private boolean isReadPairingComplete() {
		List<File> target = new ArrayList<File>();
		List<File> source = new ArrayList<File>();
		if (getContext().shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : getContext().getReference().getSequenceDictionary().getSequences()) {
				source.add(null); target.add(fsc.getAssemblyForChr(input, seq.getSequenceName()));
				source.add(null); target.add(fsc.getAssemblyMateForChr(input, seq.getSequenceName()));
			}
		} else {
			source.add(null); target.add(fsc.getAssembly(input));
			source.add(null); target.add(fsc.getAssemblyMate(input));
		}
		return IntermediateFileUtil.checkIntermediate(target, source, getContext().getConfig().ignoreFileTimestamps);
	}
	protected void process(ExecutorService threadpool) {
		if (isProcessingComplete()) return;
		log.info("START evidence assembly ", input);
		for (SAMEvidenceSource s : source) {
			if (!s.isComplete(ProcessStep.EXTRACT_READ_PAIRS) ||
				!s.isComplete(ProcessStep.EXTRACT_SOFT_CLIPS)) {
				throw new IllegalStateException(String.format("Unable to perform assembly: evidence extraction not complete for %s", s.getSourceFile()));
			}
		}
		final SAMSequenceDictionary dict = getContext().getReference().getSequenceDictionary();
		if (getContext().shouldProcessPerChromosome()) {
			final List<Callable<Void>> workers = Lists.newArrayList();
			for (int i = 0; i < dict.size(); i++) {
				final String seq = dict.getSequence(i).getSequenceName();
				
				if (IntermediateFileUtil.checkIntermediate(fsc.getAssemblyRawBamForChr(input, seq), getContext().getConfig().ignoreFileTimestamps)
					&& IntermediateFileUtil.checkIntermediate(fsc.getRealignmentFastqForChr(input, seq, 0), getContext().getConfig().ignoreFileTimestamps)) {
					log.debug("Skipping assembly for ", seq, " as output already exists.");
				} else {
					workers.add(new Callable<Void>() {
						@Override
						public Void call() {
							log.info("Starting ", seq, " breakend assembly");
							CloseableIterator<DirectedEvidence> merged = null;
							List<CloseableIterator<DirectedEvidence>> toMerge = Lists.newArrayList();
							try {
								for (SAMEvidenceSource bam : source) {
									CloseableIterator<DirectedEvidence> it = bam.iterator(
											getContext().getAssemblyParameters().includeAnomalousPairs,
											getContext().getAssemblyParameters().includeSoftClips,
											getContext().getAssemblyParameters().includeRemoteSplitReads,
											seq);
									toMerge.add(it);
								}
								merged = new AutoClosingMergedIterator<DirectedEvidence>(toMerge, DirectedEvidenceOrder.ByNatural);
								new ContigAssembler(merged, getContext().getFileSystemContext().getAssemblyRawBamForChr(input, seq), getContext().getFileSystemContext().getRealignmentFastqForChr(input, seq, 0)).run();
								merged.close();
							} catch (Exception e) {
								log.error(e, "Error performing ", seq, " breakend assembly");
								throw new RuntimeException(e);
							} finally {
								CloserUtil.close(merged);
								CloserUtil.close(toMerge);
							}
							log.info("Completed ", seq, " breakend assembly");
							return null;
						}
					});
				}
			}
			if (threadpool != null) {
				log.info("Performing multithreaded assembly");
				try {
					for (Future<Void> future : threadpool.invokeAll(workers)) {
						// Throws any exceptions back up the call stack
						try {
							future.get();
						} catch (ExecutionException e) {
							throw new RuntimeException(e);
						}
					}
				} catch (InterruptedException e) {
					String msg = "Interrupted while assembly in progress";
					log.error(e, msg);
					throw new RuntimeException(msg, e);
				}
			} else {
				log.info("Performing singlethreaded assembly");
				for (Callable<Void> c : workers) {
					try {
						c.call();
					} catch (Exception e) {
						throw new RuntimeException(e);
					}
				}
			}
		} else {
			List<CloseableIterator<DirectedEvidence>> toMerge = Lists.newArrayList();
			CloseableIterator<DirectedEvidence> merged = null;
			try {
				for (SAMEvidenceSource bam : source) {
					CloseableIterator<DirectedEvidence> it = bam.iterator(
							getContext().getAssemblyParameters().includeAnomalousPairs,
							getContext().getAssemblyParameters().includeSoftClips,
							getContext().getAssemblyParameters().includeRemoteSplitReads);
					toMerge.add(it);
				}
				merged = new AutoClosingMergedIterator<DirectedEvidence>(toMerge, DirectedEvidenceOrder.ByNatural);
				new ContigAssembler(merged, getContext().getFileSystemContext().getAssemblyRawBam(input), getContext().getFileSystemContext().getRealignmentFastq(input, 0)).run();
				merged.close();
			} finally {
				CloserUtil.close(merged);
				CloserUtil.close(toMerge);
			}
		}
		log.info("SUCCESS evidence assembly ", input);
	}
	private class ContigAssembler implements Runnable {
		private Iterator<DirectedEvidence> it;
		private File breakendOutput;
		private File realignmentFastq;
		private FastqBreakpointWriter fastqWriter = null;
		private IntervalBed throttled;
		public ContigAssembler(Iterator<DirectedEvidence> it, File breakendOutput, File realignmentFastq) {
			AssemblyConfiguration ap = getContext().getAssemblyParameters();
			this.throttled = new IntervalBed(getContext().getDictionary(), getContext().getLinear());
			DirectedEvidenceDensityThrottlingIterator dit = new DirectedEvidenceDensityThrottlingIterator(
					throttled,
					getContext().getDictionary(),
					getContext().getLinear(),
					it,
					Math.max(ap.downsampling.minimumDensityWindowSize, getMaxConcordantFragmentSize()),
					ap.downsampling.acceptDensityPortion * ap.downsampling.targetEvidenceDensity,
					ap.downsampling.targetEvidenceDensity);
			getContext().registerBuffer(breakendOutput.getName(), dit);
			this.it = dit;
			this.breakendOutput = breakendOutput;
			this.realignmentFastq = realignmentFastq;
		}
		private int generateAssembles(File unsorted) {
			SAMFileWriter writer = null;
			// Maximum value that the assembly SAM position differs from the assembly location
			int maxBreakendOffset = 0;
			try {
				long nextProgress = System.currentTimeMillis();
				SAMFileHeader samHeader = new SAMFileHeader();
				samHeader.setSequenceDictionary(getContext().getReference().getSequenceDictionary());
				samHeader.setSortOrder(SortOrder.unsorted);
				writer = getContext().getSamFileWriterFactory(false).makeSAMOrBAMWriter(samHeader, true, unsorted);
				Iterator<SAMRecordAssemblyEvidence> assemblyIt = getAssembler(it);
				SAMRecordAssemblyEvidence ass = null;
				while (assemblyIt.hasNext()) {
					ass = assemblyIt.next();
					ass = fixAssembly(ass);
					if (ass != null) {
						SAMRecord record = ass.getBackingRecord();
						maxBreakendOffset = Math.max(maxBreakendOffset, Math.abs(record.getAlignmentStart() - ass.getBreakendSummary().start));
						if (ass.getBreakendSummary() instanceof BreakpointSummary) {
							maxBreakendOffset = Math.max(maxBreakendOffset, Math.abs(record.getAlignmentStart() - ((BreakpointSummary)ass.getBreakendSummary()).start2)); 
						}
						writer.addAlignment(record);
						if (System.currentTimeMillis() > nextProgress) {
							nextProgress = System.currentTimeMillis() + 1000 * PROGRESS_UPDATE_INTERVAL_SECONDS;
							String seqname = getContext().getDictionary().getSequence(ass.getBreakendSummary().referenceIndex).getSequenceName();
							log.info(String.format("Assembly progress at %s:%d", seqname, ass.getBreakendSummary().start));
							if (getContext().getConfig().getVisualisation().assembly) {
								File exportDir = getContext().getConfig().getVisualisation().directory;
								if (assemblyIt instanceof PositionalAssembler) {
									if (exportDir != null) {
										File exportFilename = new File(exportDir, String.format("assembly_%s_%d.dot", seqname, ass.getBreakendSummary().start)); 
										((PositionalAssembler)assemblyIt).exportGraph(exportFilename);
									}
								} else if (assemblyIt instanceof DirectedPositionalAssembler) {
									if (exportDir != null) {
										((DirectedPositionalAssembler)assemblyIt).exportGraph(
												new File(exportDir, String.format("assemblyf_%s_%d.dot", seqname, ass.getBreakendSummary().start)),
												new File(exportDir, String.format("assemblyb_%s_%d.dot", seqname, ass.getBreakendSummary().start)));
									}
								}
							}
						}
					}
				}
				if (ass != null) {
					log.info(String.format("Assembly complete for %s", getContext().getDictionary().getSequence(ass.getBreakendSummary().referenceIndex).getSequenceName()));
				}
				try {
					writer.close();
				} finally {
					writer = null;
				}
				String throttleFilename = breakendOutput.getAbsolutePath() + ".throttled.bed";
				try {
					throttled.write(new File(throttleFilename), "Regions of high coverage where only a portion of supporting reads considered for assembly");
				} catch (IOException e) {
					log.warn(e, "Unable to write " + throttleFilename);
				}
			} finally {
				CloserUtil.close(writer);
			}
			return maxBreakendOffset;
		}
		private SAMRecordAssemblyEvidence fixAssembly(SAMRecordAssemblyEvidence ass) {
			if (ass == null) return null;
			AssemblyConfiguration ap = getContext().getAssemblyParameters();
			// realign
			byte[] be = ass.getBreakendSequence();
			if (be != null && be.length > ap.maxExpectedBreakendLengthMultiple * getMaxConcordantFragmentSize()) {
				log.debug(String.format("Filtering %s at %s due to misassembly (breakend %dbp)",
						ass.getEvidenceID(),
						ass.getBreakendSummary(),
						be.length));
				return null;
			}
    		if (ap.anchorRealignment.perform && ass.isBreakendExact()) {
    			int realignmentWindowSize = (int)(ap.anchorRealignment.realignmentWindowReadLengthMultiples * getMaxReadLength());
    			SAMRecordAssemblyEvidence fullRealignment = ass.realign(realignmentWindowSize, true, ap.anchorRealignment.realignmentMinimumAnchorRetainment);
    			// use full assembly realignment to find small indels and reference assemblies
    			if (fullRealignment != null && (fullRealignment.isReferenceAssembly() || fullRealignment.isSpanningAssembly())) {
    				ass = fullRealignment;
    			} else {
    				// use targeted realignment to find the correct breakend location
    				ass = ass.realign(realignmentWindowSize, false, ap.anchorRealignment.realignmentMinimumAnchorRetainment);
    			}
    		}
			if (ass == null || ass.getBreakendSummary() == null) return null;
			getContext().getAssemblyParameters().applyBasicFilters(ass);
			if (getContext().getAssemblyParameters().writeFiltered) return ass;
			if (ass.isReferenceAssembly()) return null;
			if (ass.isAssemblyFiltered()) return null;
			return ass;
		}
		private void assembliesToSortedBamAndFastq(int maxBreakendOffset, File unsorted) throws IOException {
			File tmp = FileSystemContext.getWorkingFileFor(breakendOutput);
			File fq = FileSystemContext.getWorkingFileFor(realignmentFastq);
			fq.delete();
			tmp.delete();
			CloseableIterator<SAMRecord> readerIt = null;
			try {
				// sort the assemblies
				new SortCallable(getContext(), unsorted, breakendOutput, SortOrder.coordinate, header -> {
					header.addComment(String.format("gridss.maxBreakendOffset=%d", maxBreakendOffset));
					return header;
				}).call();
				log.info("Writing assembly fastq " + fq);
				// then write the fastq
				fastqWriter = new FastqBreakpointWriter(getContext().getFastqWriterFactory().newWriter(fq));
				readerIt = getContext().getSamReaderIterator(breakendOutput, SortOrder.coordinate);
				while (readerIt.hasNext()) {
					SAMRecord record = readerIt.next();
					SAMRecordAssemblyEvidence assembly = AssemblyFactory.hydrate(AssemblyEvidenceSource.this, record);
					if (getContext().getRealignmentParameters().shouldRealignBreakend(assembly)) {
						fastqWriter.write(assembly);
					}
				}
				readerIt.close();
				readerIt = null;
				fastqWriter.close();
				fastqWriter = null;
				FileHelper.move(fq, realignmentFastq, true);
			} finally {
				CloserUtil.close(readerIt);
				CloserUtil.close(fastqWriter);
			}
		}
		@Override
		public void run() {
			try {
				File unsorted = FileSystemContext.getWorkingFileFor(breakendOutput, "unsorted."); 
				unsorted.delete();
				int maxBreakendOffset = generateAssembles(unsorted);
				assembliesToSortedBamAndFastq(maxBreakendOffset, unsorted);
			} catch (Throwable e) {
				log.error(e, "Error assembling breakend ", breakendOutput);
				e.printStackTrace();
				throw new RuntimeException("Error assembling breakend", e);
			}
		}
	}
	protected Iterator<SAMRecordAssemblyEvidence> getAssembler(Iterator<DirectedEvidence> it) {
		switch (getContext().getAssemblyParameters().method) {
			case Positional:
				return new DirectedPositionalAssembler(getContext(), this, it);
				//return new PositionalAssembler(getContext(), this, it);
			case Subgraph:
				return new ReadEvidenceAssemblyIterator(new DeBruijnSubgraphAssembler(getContext(), this), it);
		}
		throw new IllegalArgumentException("Assembly algorithm has not been set");
    }
	public int getAssemblyEvidenceWindowSize() {
		return (int)(getContext().getAssemblyParameters().subgraph.assemblyMargin * getMaxConcordantFragmentSize());
	}
	public int getAssemblyMaximumEvidenceDelay() {
		return Math.max(getContext().getAssemblyParameters().subgraph.minSubgraphWidthForTimeout,
				(int)(getContext().getAssemblyParameters().subgraph.maxSubgraphFragmentWidth * getMaxConcordantFragmentSize()));
	}
	public int getAssemblyWindowSize() {
		// Positional assembly should have a smaller window size
		// but the width of the assembly itself is unbounded
		// in the degenerate misassembly case
		return getAssemblyMaximumEvidenceDelay() + 3 * getAssemblyEvidenceWindowSize() + 2;
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
}
