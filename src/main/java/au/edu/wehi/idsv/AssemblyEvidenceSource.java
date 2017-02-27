package au.edu.wehi.idsv;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import com.google.common.base.Stopwatch;
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.util.concurrent.MoreExecutors;

import au.edu.wehi.idsv.SequentialEvidenceAssemblyAllocator.BreakendAssemblyEvidenceSupport;
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
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.sam.GatherBamFiles;

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
	private final IntervalBed throttled;
	/**
	 * Generates assembly evidence based on the given evidence
	 * @param evidence evidence for creating assembly
	 * @param intermediateFileLocation location to store intermediate files
	 */
	public AssemblyEvidenceSource(ProcessingContext processContext, List<SAMEvidenceSource> evidence, File assemblyFile) {
		super(processContext, assemblyFile, null, -1);
		this.source = evidence;
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
		List<QueryInterval> chunks = getContext().getReference().getIntervals(getContext().getConfig().chunkSize);
		List<File> assembledChunk = new ArrayList<>();
		List<Future<Void>> tasks = new ArrayList<>();
		for (int i = 0; i < chunks.size(); i++) {
			QueryInterval chunck = chunks.get(i);
			File f = getContext().getFileSystemContext().getAssemblyChunkBam(getFile(), i);
			int chunkNumber = i;
			assembledChunk.add(f);
			if (!f.exists()) {
				tasks.add(threadpool.submit(() -> { assembleChunk(f, chunkNumber, chunck); return null; }));
			}
			
		}
		runTasks(tasks);
		log.info("Breakend assembly complete.");
		List<File> deduplicatedChunks = assembledChunk;
		if (Iterables.any(source, ses -> ses.getMetrics().getIdsvMetrics().SECONDARY_NOT_SPLIT > 0)) {
			if (!getContext().getConfig().multimapping) {
				log.info("Not performing unique assembly allocation of multimapping reads as the configuration setting multimapping is set to false.");
			} else if (!getContext().getConfig().multimappingUniqueAssemblyAllocation) {
				log.info("Not performing unique assembly allocation of multimapping reads as the configuration setting multimappingUniqueAssemblyAllocation is set to false.");
			} else { 
				log.info("Multimapping mode invoked due to existiance of at least one BAM file with a non-split secondary alignment.");
				log.info("Uniquely assigning multi-mapping reads to alignment location supported by the best assembly.");
				log.info("Loading assembly evidence allocations");
				deduplicatedChunks = assembledChunk.stream()
						.map(f -> FileSystemContext.getWorkingFileFor(f, "deduplicated."))
						.collect(Collectors.toList());
				// technically this is an overestimation as a single alignment can
				// a DP or OEA and an Indel/SC/SR.
				long svReadAlignmentCount = source.stream()
					.map(ses -> ses.getSVMetrics())
					.mapToLong(m -> m.DISCORDANT_READ_PAIR_ALIGNMENTS +
							m.UNMAPPED_MATE_READ_ALIGNMENTS +
							m.STRUCTURAL_VARIANT_READ_ALIGNMENTS)
					.sum();
				try (GreedyAssemblyAllocationCache cache = new GreedyAssemblyAllocationCache(svReadAlignmentCount)) {
					tasks = new ArrayList<>();
					for (int i = 0; i < chunks.size(); i++) {
						QueryInterval chunk = chunks.get(i);
						File in = assembledChunk.get(i);
						tasks.add(threadpool.submit(() -> { loadAssemblyEvidenceAllocation(cache, in, chunk); return null; }));
					}
					runTasks(tasks);
					
					log.info("Allocating multi-mapping reads to assemblies");
					tasks = new ArrayList<>();
					for (int i = 0; i < chunks.size(); i++) {
						QueryInterval chunk = chunks.get(i);
						File in = assembledChunk.get(i);
						File out = deduplicatedChunks.get(i);
						if (!out.exists()) {
							tasks.add(threadpool.submit(() -> { deduplicateChunk(in, out, chunk, cache); return null; }));
						}
					}
					runTasks(tasks);
				}
			}
		}
		log.info("Merging assembly files");
		// Merge chunk files
		File tmpout = FileSystemContext.getWorkingFileFor(getFile());
		GatherBamFiles gather = new picard.sam.GatherBamFiles();
		List<String> args = new ArrayList<>();
		for (File f : deduplicatedChunks) {
			args.add("INPUT=" + f.getAbsolutePath());
		}
		args.add("OUTPUT=" + tmpout.getAbsolutePath());
		int returnCode = gather.instanceMain(args.toArray(new String[] {}));
		if (returnCode != 0) {
			String msg = String.format("Error executing GatherBamFiles. GatherBamFiles returned status code %d", returnCode);
			log.error(msg);
			throw new RuntimeException(msg);
		}
		// Sorting is not required since each chunk was already sorted, and each chunk
		// contains sequential genomic coordinates. We also don't need to index as we only need assembly.sv.bam indexed
		// SAMFileUtil.sort(getContext().getFileSystemContext(), tmpout, getFile(), SortOrder.coordinate);
		FileHelper.move(tmpout, getFile(), true);
		if (gridss.Defaults.DELETE_TEMPORARY_FILES) {
			FileHelper.delete(tmpout, true);
			for (File f : assembledChunk) {
				FileHelper.delete(f, true);
			}
			for (File f : deduplicatedChunks) {
				FileHelper.delete(f, true);
			}
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
	private void runTasks(List<Future<Void>> tasks) {
		// Assemble as much as we can before dying
		Exception firstException = null;
		for (Future<Void> f : tasks) {
			try {
				f.get();
			} catch (Exception e) {
				if (firstException == null) {
					firstException = e;
				}
			}
		}
		log.info("Breakend assembly complete.");
		if (firstException != null) {
			log.error(firstException, "Fatal error during assembly ");
			throw new RuntimeException(firstException);
		}
	}
	private void assembleChunk(File output, int chunkNumber, QueryInterval qi) throws IOException {
		AssemblyIdGenerator assemblyNameGenerator = new SequentialIdGenerator(String.format("asm%d-", chunkNumber));
		String chuckName = String.format("chunk %d (%s:%d-%d)", chunkNumber, getContext().getReference().getSequenceDictionary().getSequence(qi.referenceIndex).getSequenceName(), qi.start, qi.end);
		log.info(String.format("Starting assembly on %s", chuckName));
		Stopwatch timer = Stopwatch.createStarted();
		SAMFileHeader header = getContext().getBasicSamHeader();
		// TODO: add assembly @PG header
		File filteredout = FileSystemContext.getWorkingFileFor(output, "filtered.");
		File tmpout = FileSystemContext.getWorkingFileFor(output, "gridss.tmp.");
		try (SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, tmpout)) {
			if (getContext().getAssemblyParameters().writeFiltered) {
				try (SAMFileWriter filteredWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, filteredout)) {
					for (BreakendDirection direction : BreakendDirection.values()) {
						assembleChunk(writer, filteredWriter, qi, direction, assemblyNameGenerator);
					}
				}
			} else {
				for (BreakendDirection direction : BreakendDirection.values()) {
					assembleChunk(writer, null, qi, direction, assemblyNameGenerator);
				}
			}
		} catch (Exception e) {
			log.error(e, "Error assembling ", chuckName);
			if (getContext().getConfig().terminateOnFirstError) {
				System.exit(1);
			}
			throw e;
		} finally {
			timer.stop();
			log.info(String.format("Completed assembly on %s in %ds (%s)", chuckName, timer.elapsed(TimeUnit.SECONDS), timer.toString()));
		}
		SAMFileUtil.sort(getContext().getFileSystemContext(), tmpout, output, SortOrder.coordinate);
		if (gridss.Defaults.DELETE_TEMPORARY_FILES) {
			tmpout.delete();
			filteredout.delete();
		}
	}
	private QueryInterval getExpanded(QueryInterval qi) {
		int expansion = (int)(2 * getMaxConcordantFragmentSize() * getContext().getConfig().getAssembly().maxExpectedBreakendLengthMultiple) + 1;
		int end = getContext().getReference().getSequenceDictionary().getSequence(qi.referenceIndex).getSequenceLength();
		QueryInterval expanded = new QueryInterval(qi.referenceIndex, Math.max(1, qi.start - expansion), Math.min(end, qi.end + expansion));
		return expanded;
	}
	private void assembleChunk(SAMFileWriter writer, SAMFileWriter filteredWriter, QueryInterval qi, BreakendDirection direction, AssemblyIdGenerator assemblyNameGenerator) {
		QueryInterval expanded = getExpanded(qi);
		try (CloseableIterator<DirectedEvidence> input = mergedIterator(source, expanded)) {
			Iterator<DirectedEvidence> throttledIt = throttled(input);
			Iterator<SAMRecord> evidenceIt = new PositionalAssembler(getContext(), AssemblyEvidenceSource.this, assemblyNameGenerator, throttledIt, direction);
			while (evidenceIt.hasNext()) {
				SAMRecord asm = evidenceIt.next();
				asm = transformAssembly(asm); // transform before chunk bounds checking as the position may have moved
				if (asm.getAlignmentStart() >= qi.start && asm.getAlignmentStart() <= qi.end) {
					// only output assemblies that start within our chunk
					if (shouldFilterAssembly(asm)) {
						if (filteredWriter != null) {
							filteredWriter.addAlignment(asm);
						}
					} else {
						writer.addAlignment(asm);
					}
				}
			}
		}
	}
	private void loadAssemblyEvidenceAllocation(GreedyAssemblyAllocationCache cache, File in, QueryInterval qi) throws IOException {
		log.debug(String.format("Caching assembly evidence allocation in interval %s:%d-%d", getContext().getDictionary().getSequence(qi.referenceIndex).getSequenceName(), qi.start, qi.end));
		QueryInterval expanded = getExpanded(qi);
		try (CloseableIterator<DirectedEvidence> reads = mergedIterator(source, expanded)) {
			try (SamReader reader = factory.open(in)) {
				try (SAMRecordIterator assemblies = reader.iterator()) {
					SequentialEvidenceAssemblyAllocator it = new SequentialEvidenceAssemblyAllocator(
							getContext().getLinear(),
							reads,
							assemblies, 0);
					while (it.hasNext()) {
						BreakendAssemblyEvidenceSupport allocation = it.next();
						float score = 0;
						for (DirectedEvidence assemblyEvidence : SingleReadEvidence.createEvidence(this, 0, allocation.assemblyRecord)) {
							score = Math.max(score, assemblyEvidence.getBreakendQual());
						}
						for (DirectedEvidence support : allocation.support) {
							if (support.isFromMultimappingFragment()) {
								cache.addBreakendAssemblyAllocation(score, support);
							}
						}
					}
				}
			}
		}
	}
	private void deduplicateChunk(File in, File out, QueryInterval qi, GreedyAssemblyAllocationCache cache) throws IOException {
		log.debug(String.format("Uniquely assigning assembly evidence allocation in interval %s:%d-%d", getContext().getDictionary().getSequence(qi.referenceIndex).getSequenceName(), qi.start, qi.end));
		File tmpout = FileSystemContext.getWorkingFileFor(out);
		QueryInterval expanded = getExpanded(qi);
		try (CloseableIterator<DirectedEvidence> reads = mergedIterator(source, expanded)) {
			try (SamReader reader = factory.open(in)) {
				try (SAMRecordIterator assemblies = reader.iterator()) {
					try (SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), true, tmpout)) {
						SequentialEvidenceAssemblyAllocator it = new SequentialEvidenceAssemblyAllocator(
								getContext().getLinear(),
								reads,
								assemblies, 0);
						while (it.hasNext()) {
							BreakendAssemblyEvidenceSupport allocation = it.next();
							SAMRecord asm = allocation.assemblyRecord;
							List<DirectedEvidence> newSupport = new ArrayList<>(allocation.support.size());
							for (DirectedEvidence support : allocation.support) {
								if (!support.isFromMultimappingFragment() || cache.isBestBreakendAssemblyAllocation(support)) {
									newSupport.add(support);
								}
							}
							AssemblyAttributes.annotateAssembly(getContext(), asm, newSupport);
							if (shouldFilterAssembly(asm)) {
								// TODO: append to filtered bam
							} else {
								writer.addAlignment(asm);
							}
						}
					}
				}
			}
		}
		FileHelper.move(tmpout, out, true);
	}
	@Override
	public void ensureExtracted() throws IOException {
		ensureMetrics();
		File svFile = getContext().getFileSystemContext().getSVBam(getFile());
		File withsplitreadsFile = FileSystemContext.getWorkingFileFor(svFile, "gridss.tmp.withsplitreads.");
		ensureMetrics();
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
		List<SingleReadEvidence> breakends = SingleReadEvidence.createEvidence(this, 0, asm);
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
	protected int getSortWindowSize() {
		// Worse case for assembly is worse than for the input BAMs as
		// assembly allows for unanchored split reads using the XNX notation
		int windowSize = super.getSortWindowSize() + getMaxConcordantFragmentSize();
		return windowSize;
	}
	@Override
	public int getMaxConcordantFragmentSize() {
		int maxSourceFragSize = source.stream().mapToInt(s -> s.getMaxConcordantFragmentSize()).max().orElse(0);
		if (getFile().exists()) {
			maxSourceFragSize = Math.max(super.getMaxConcordantFragmentSize(), maxSourceFragSize);
		}
		return maxSourceFragSize;
	}
	@Override
	public int getMinConcordantFragmentSize() {
		int minSourceFragSize = source.stream().mapToInt(s -> s.getMinConcordantFragmentSize()).max().orElse(0);
		if (getFile().exists()) {
			minSourceFragSize = Math.max(super.getMaxConcordantFragmentSize(), minSourceFragSize);
		}
		return minSourceFragSize;
	}
	/**
	 * Maximum read length of reads contributing to assemblies
	 * @return
	 */
	@Override
	public int getMaxReadLength() {
		int maxReadLength = source.stream().mapToInt(s -> s.getMaxReadLength()).max().orElse(0);
		if (getFile().exists()) {
			maxReadLength = Math.max(super.getMaxReadLength(), maxReadLength);
		}
		return maxReadLength;
	}
	@Override
	public int getMaxReadMappedLength() {
		int maxMappedReadLength = source.stream().mapToInt(s -> s.getMaxReadMappedLength()).max().orElse(0);
		if (getFile().exists()) {
			maxMappedReadLength = Math.max(super.getMaxReadMappedLength(), maxMappedReadLength);
		}
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
			Iterator<SAMRecord> evidenceIt = new PositionalAssembler(getContext(), this, new SequentialIdGenerator("asm"), loggedIt, direction);
	    	list.add(evidenceIt);
		}
		return Iterators.concat(list.iterator());
	}
	/**
	 * Assembly contigs are not paired
	 */
	@Override
	public boolean knownSingleEnded() { return true; } 
}
