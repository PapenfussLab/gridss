package au.edu.wehi.idsv;

import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.configuration.AssemblyConfiguration;
import au.edu.wehi.idsv.debruijn.positional.PositionalAssembler;
import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.sam.SAMFileUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.sam.SamTags;
import au.edu.wehi.idsv.util.FileHelper;
import au.edu.wehi.idsv.visualisation.AssemblyTelemetry;
import com.google.common.base.Stopwatch;
import com.google.common.util.concurrent.MoreExecutors;
import gridss.SoftClipsToSplitReads;
import gridss.cmdline.CommandLineProgramHelper;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

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
	public static final String INPUT_CATEGORY_SAM_HEADER_PREFIX = "gridss_input_category=";
	private final List<SAMEvidenceSource> source;
	private int cachedMaxSourceFragSize = -1;
	private int cachedMinConcordantFragmentSize = -1;
	private int cachedMaxReadLength = -1;
	private int cachedMaxReadMappedLength = -1;
	private AssemblyTelemetry telemetry;
	private SAMFileHeader header;
	protected List<String> assembledCategories;
	private int[] assemblyOrdinalToProcessingCategoryLookup;
	/**
	 * Generates assembly evidence based on the given evidence
	 * @param evidence evidence for creating assembly
	 */
	public AssemblyEvidenceSource(ProcessingContext processContext, List<SAMEvidenceSource> evidence, File assemblyFile) {
		super(processContext, assemblyFile, null, -1);
		this.source = filterOutLongReadLibraries(evidence);
	}

	private static List<SAMEvidenceSource> filterOutLongReadLibraries(List<SAMEvidenceSource> evidence) {
		List<SAMEvidenceSource> filtered = evidence.stream().
				filter(e -> !e.isLongReadLibrary())
				.collect(Collectors.toList());
		return filtered;
	}

	public void assembleBreakends(ExecutorService threadpool) throws IOException {
		assembleBreakends(threadpool, 0, 1);
	}
	/**
	 * Perform breakend assembly 
	 * @param threadpool 
	 * @throws IOException 
	 */
	public void assembleBreakends(ExecutorService threadpool, int jobNodeIndex, int jobNodes) throws IOException {
		IntervalBed excludedRegions = new IntervalBed(getContext().getLinear());
		IntervalBed safetyRegions = new IntervalBed(getContext().getLinear());
		IntervalBed downsampledRegions = new IntervalBed(getContext().getLinear());
		source.stream().forEach(ses -> ses.assertPreprocessingComplete());
		if (threadpool == null) {
			threadpool = MoreExecutors.newDirectExecutorService();
		}
		invalidateSummaryCache();
		if (getContext().getConfig().getVisualisation().assemblyTelemetry) {
			telemetry = new AssemblyTelemetry(getContext().getFileSystemContext().getAssemblyTelemetry(getFile(), jobNodeIndex), getContext().getDictionary());
		}
		List<QueryInterval[]> chunks = getContext().getReference().getIntervals(getContext().getConfig().chunkSize, getContext().getConfig().chunkSequenceChangePenalty);
		List<File> assembledChunk = new ArrayList<>();
		List<Future<Void>> tasks = new ArrayList<>();
		for (int i = 0; i < chunks.size(); i++) {
			if (i % jobNodes == jobNodeIndex) {
				QueryInterval[] chunk = chunks.get(i);
				File f = getContext().getFileSystemContext().getAssemblyChunkBam(getFile(), i);
				int chunkNumber = i;
				assembledChunk.add(f);
				if (!f.exists()) {
					tasks.add(threadpool.submit(() -> {
						assembleChunk(f, chunkNumber, chunk, excludedRegions, safetyRegions, downsampledRegions);
						return null;
					}));
				}
			}
		}
		runTasks(tasks);
		if (telemetry != null) {
			telemetry.close();
			telemetry = null;
		}
		excludedRegions.write(getContext().getFileSystemContext().getAssemblyExcludedRegions(getFile(), jobNodeIndex), "excludedDueToGraphComplexity");
		safetyRegions.write(getContext().getFileSystemContext().getAssemblySafetyRegions(getFile(), jobNodeIndex), "subsetOfContigsCalledDueToGraphComplexity");
		downsampledRegions.write(getContext().getFileSystemContext().getAssemblyDownsampledRegions(getFile(), jobNodeIndex), "subsetOfReadsAssembled");
		log.info(String.format("Breakend assembly complete (node %d, %d total)", jobNodeIndex, jobNodes));
		List<File> deduplicatedChunks = assembledChunk;
		long secondaryNotSplit = source.stream().mapToLong(ses -> ses.getMetrics().getIdsvMetrics().SECONDARY_NOT_SPLIT).sum();
		if (secondaryNotSplit > 0) {
			log.warn(String.format("Found %d secondary alignments that were not split read alignments. GRIDSS no longer supports multi-mapping alignment. These reads will be ignored.", secondaryNotSplit));
		}
		if (jobNodes > 1) {
			log.info("Not merging assembly files since not all chunks were assembled.");
			return;
		}
		log.info("Merging assembly files");
		// Merge chunk files
		File out = getFile();
		File tmpout = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(getFile()) : out;
		CommandLineProgramHelper gather = new CommandLineProgramHelper(new picard.sam.GatherBamFiles());
		for (File f : deduplicatedChunks) {
			gather.addArg("INPUT", f.getPath());
		}
		gather.addArg("OUTPUT", tmpout.getPath());
		if (getContext().getCommandLineProgram() != null) {
			gather.setCommonArgs(getContext().getCommandLineProgram());
		}
		int returnCode = gather.run();
		if (returnCode != 0) {
			String msg = String.format("Error executing GatherBamFiles. GatherBamFiles returned status code %d", returnCode);
			log.error(msg);
			throw new RuntimeException(msg);
		}
		// Sorting is not required since each chunk was already sorted, and each chunk
		// contains sequential genomic coordinates. We also don't need to index as we only need assembly.sv.bam indexed
		// SAMFileUtil.sort(getContext().getFileSystemContext(), tmpout, getFile(), SortOrder.coordinate);
		if (tmpout != out) {
			FileHelper.move(tmpout, out, true);
		}
		invalidateSummaryCache();
		if (gridss.Defaults.DELETE_TEMPORARY_FILES) {
			if (tmpout != out) {
				FileHelper.delete(tmpout, true);
			}
			for (File f : assembledChunk) {
				FileHelper.delete(f, true);
			}
			for (File f : deduplicatedChunks) {
				FileHelper.delete(f, true);
			}
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
		if (firstException != null) {
			log.error(firstException, "Fatal error during assembly ");
			throw new RuntimeException(firstException);
		}
		log.info("Breakend assembly complete.");
	}
	private void assembleChunk(File output, int chunkNumber, QueryInterval[] qi, IntervalBed excludedRegions, IntervalBed safetyRegions, IntervalBed downsampledRegions) throws IOException {
		AssemblyIdGenerator assemblyNameGenerator = new SequentialIdGenerator(String.format(getContext().getConfig().getAssembly().contigNamePrefix, chunkNumber));
		String chuckName = String.format("chunk %d (%s:%d-%s:%d)", chunkNumber,
			getContext().getDictionary().getSequence(qi[0].referenceIndex).getSequenceName(), qi[0].start,
			getContext().getDictionary().getSequence(qi[qi.length-1].referenceIndex).getSequenceName(), qi[qi.length-1].end);
		log.info(String.format("Starting assembly on %s", chuckName));
		Stopwatch timer = Stopwatch.createStarted();
		File filteredout = FileSystemContext.getWorkingFileFor(output, "filtered.");
		File tmpout = FileSystemContext.getWorkingFileFor(output, "gridss.tmp.");
		try (SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(getHeader(), false, tmpout)) {
			if (getContext().getAssemblyParameters().writeFiltered) {
				try (SAMFileWriter filteredWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(getHeader(), false, filteredout)) {
					for (BreakendDirection direction : BreakendDirection.values()) {
						assembleChunk(writer, filteredWriter, chunkNumber, qi, direction, assemblyNameGenerator, excludedRegions, safetyRegions, downsampledRegions);
					}
				}
			} else {
				for (BreakendDirection direction : BreakendDirection.values()) {
					assembleChunk(writer, null, chunkNumber, qi, direction, assemblyNameGenerator, excludedRegions, safetyRegions, downsampledRegions);
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
		if (gridss.Defaults.DEFENSIVE_GC) {
			log.info("Requesting defensive GC to ensure OS file handles are closed");
			System.gc();
			System.runFinalization();
		}
	}

	private QueryInterval[] getExpanded(QueryInterval[] intervals) {
		QueryInterval[] expanded = QueryIntervalUtil.padIntervals(
				getContext().getDictionary(),
				intervals,
				// expand bounds to keep any contig that could overlap our intervals
				(int)(2 * getMaxConcordantFragmentSize() * getContext().getConfig().getAssembly().maxExpectedBreakendLengthMultiple) + 1);
		return expanded;
	}
	private void assembleChunk(SAMFileWriter writer, SAMFileWriter filteredWriter, int chunkNumber, QueryInterval[] intervals, BreakendDirection direction, AssemblyIdGenerator assemblyNameGenerator,
							   IntervalBed excludedRegions, IntervalBed safetyRegions, IntervalBed downsampledRegions) {
		QueryInterval[] expanded = getExpanded(intervals);
		try (CloseableIterator<DirectedEvidence> input = mergedIterator(source, expanded, EvidenceSortOrder.SAMRecordStartPosition)) {
			Iterator<DirectedEvidence> throttledIt = throttled(input, downsampledRegions);
			PositionalAssembler assembler = new PositionalAssembler(getContext(), AssemblyEvidenceSource.this, assemblyNameGenerator, throttledIt, direction, excludedRegions, safetyRegions);
			if (telemetry != null) {
				assembler.setTelemetry(telemetry.getTelemetry(chunkNumber, direction));
			}
			while (assembler.hasNext()) {
				SAMRecord asm = assembler.next();
				asm = transformAssembly(asm); // transform before chunk bounds checking as the position may have moved
				if (QueryIntervalUtil.overlaps(intervals, asm.getReferenceIndex(), asm.getAlignmentStart())) {
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
	@Override
	public synchronized void ensureExtracted() throws IOException {
		ensureMetrics();
		File svFile = getContext().getFileSystemContext().getSVBam(getFile());
		File withsplitreadsFile = FileSystemContext.getWorkingFileFor(svFile, "gridss.tmp.withsplitreads.");
		ensureMetrics();
		if (!svFile.exists()) {
			log.info("Identifying split reads for " + getFile().getAbsolutePath());
			SoftClipsToSplitReads program = new SoftClipsToSplitReads();
			program.setReference(getProcessContext().getReference());
			program.setFileSystemContext(getProcessContext().getFileSystemContext());
			CommandLineProgramHelper cmd = new CommandLineProgramHelper(program);
			cmd.addArg("WORKER_THREADS", getProcessContext().getWorkerThreadCount());
			cmd.addArg("INPUT", getFile().getPath());
			cmd.addArg("OUTPUT", svFile.getPath());
			cmd.addArg("READJUST_PRIMARY_ALIGNMENT_POSITION", "true");
			cmd.addArg("REALIGN_ENTIRE_READ=", getContext().getConfig().getAssembly().realignContigs);
			execute(cmd);
		}
		SAMFileUtil.sort(getContext().getFileSystemContext(), withsplitreadsFile, svFile, SortOrder.coordinate);
	}
	@Override
	public boolean shouldFilter(SAMRecord r) {
		if (r.hasAttribute("OA") && !SAMRecordUtil.overlapsOriginalAlignment(r)) {
			return true;
		}
		int[] supportCategory = r.getSignedIntArrayAttribute(SamTags.ASSEMBLY_EVIDENCE_CATEGORY);
		if (supportCategory != null && Arrays.stream(supportCategory).anyMatch(c -> c >= getContext().getCategoryCount())) {
			String msg = String.format("Fatal error: Assembly contig %s contains for unknown input. " +
					"GRIDSS performs joint assembly and does not support per-input assembling. " +
					"inputs and labels must be the same for all assembly and variant calling steps. ",
					r.getReadName());
			log.error(msg);
			throw new RuntimeException(msg);
		}
		return super.shouldFilter(r);
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
		if (attr.getOriginatingFragmentID(null, null, null, null).size() < ap.minReads) {
			return true;
		}
		// unanchored assembly that not actually any longer than any of the reads that were assembled
		if (AssemblyAttributes.isUnanchored(asm) && breakendLength <= attr.getAssemblyMaxReadLength()) {
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
	private Iterator<DirectedEvidence> throttled(Iterator<DirectedEvidence> it, IntervalBed downsampledRegions) {
		AssemblyConfiguration ap = getContext().getAssemblyParameters();
		DirectedEvidenceDensityThrottlingIterator dit = new DirectedEvidenceDensityThrottlingIterator(
				downsampledRegions,
				getContext().getDictionary(),
				getContext().getLinear(),
				it,
				EvidenceSortOrder.SAMRecordStartPosition,
				Math.max(ap.downsampling.minimumDensityWindowSize, getMaxConcordantFragmentSize()),
				ap.downsampling.acceptDensityPortion * ap.downsampling.targetEvidenceDensity,
				ap.downsampling.targetEvidenceDensity,
				true, false);
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
	private void invalidateSummaryCache() {
		cachedMaxSourceFragSize = -1;
		cachedMinConcordantFragmentSize = -1;
		cachedMaxReadLength = -1;
		cachedMaxReadMappedLength = -1;
	}
	private int calcMaxConcordantFragmentSize() {
		int fs = source.stream().mapToInt(s -> s.getMaxConcordantFragmentSize()).max().orElse(0);
		if (getFile().exists()) {
			fs = Math.max(super.getMaxConcordantFragmentSize(), fs);
		}
		return fs;
	}
	@Override
	public int getMaxConcordantFragmentSize() {
		if (cachedMaxSourceFragSize == -1) {
			cachedMaxSourceFragSize = calcMaxConcordantFragmentSize();
		}
		return cachedMaxSourceFragSize;
	}
	@Override
	public int getMinConcordantFragmentSize() {
		if (cachedMinConcordantFragmentSize == -1) {
			cachedMinConcordantFragmentSize = calcMinConcordantFragmentSize();
		}
		return cachedMinConcordantFragmentSize;
	}
	public int calcMinConcordantFragmentSize() {
		int minSourceFragSize = source.stream().mapToInt(s -> s.getMinConcordantFragmentSize()).min().orElse(0);
		if (getFile().exists()) {
			minSourceFragSize = Math.min(super.getMinConcordantFragmentSize(), minSourceFragSize);
		}
		return minSourceFragSize;
	}
	/**
	 * Maximum read length of reads contributing to assemblies
	 * @return
	 */
	@Override
	public int getMaxReadLength() {
		if (cachedMaxReadLength == -1) {
			cachedMaxReadLength = calcMaxReadLength();
		}
		return cachedMaxReadLength;
	}
	public int calcMaxReadLength() {
		int maxReadLength = source.stream().mapToInt(s -> s.getMaxReadLength()).max().orElse(0);
		if (getFile().exists()) {
			maxReadLength = Math.max(super.getMaxReadLength(), maxReadLength);
		}
		return maxReadLength;
	}
	@Override
	public int getMaxReadMappedLength() {
		if (cachedMaxReadMappedLength == -1) {
			cachedMaxReadMappedLength = calcMaxReadMappedLength();
		}
		return cachedMaxReadMappedLength;
	}
	public int calcMaxReadMappedLength() {
		int maxMappedReadLength = source.stream().mapToInt(s -> s.getMaxReadMappedLength()).max().orElse(0);
		if (getFile().exists()) {
			maxMappedReadLength = Math.max(super.getMaxReadMappedLength(), maxMappedReadLength);
		}
		return maxMappedReadLength;
	}
	/*
	@SuppressWarnings("unused")
	private Iterator<SAMRecord> getAllAssemblies_single_threaded() {
		ProgressLogger progressLog = new ProgressLogger(log);
		List<Iterator<SAMRecord>> list = new ArrayList<>();
		for (BreakendDirection direction : BreakendDirection.values()) {
			CloseableIterator<DirectedEvidence> it = mergedIterator(source, false, EvidenceSortOrder.SAMRecordStartPosition);
			Iterator<DirectedEvidence> throttledIt = throttled(it, new IntervalBed(getContext().getLinear()));
			ProgressLoggingDirectedEvidenceIterator<DirectedEvidence> loggedIt = new ProgressLoggingDirectedEvidenceIterator<>(getContext(), throttledIt, progressLog);
			Iterator<SAMRecord> evidenceIt = new PositionalAssembler(getContext(), this, new SequentialIdGenerator("asm"), loggedIt, direction, new IntervalBed(getContext().getLinear()), new IntervalBed(getContext().getLinear()));
	    	list.add(evidenceIt);
		}
		return Iterators.concat(list.iterator());
	}
	*/
	/**
	 * Assembly contigs are not paired
	 */
	@Override
	public boolean knownSingleEnded() { return true; }

	@Override
	protected SamReader getReader() {
		SamReader reader = super.getReader();
		this.header = reader.getFileHeader();
		this.assembledCategories = getAssemblyCategories(getHeader());
		return reader;
	}

	/**
	 * Checks that all input categories have been assembled.
	 * @param aesList
	 */
	public static void validateAllCategoriesAssembled(ProcessingContext pc, List<AssemblyEvidenceSource> aesList) {
		Set<String> pcCategories = IntStream.range(0, pc.getCategoryCount())
				.mapToObj(x -> pc.getCategoryLabel(x))
				.collect(Collectors.toSet());
		for (int i = 0; i < aesList.size(); i++) {
			AssemblyEvidenceSource aes1 = aesList.get(i);
			// Assembled but not in input
			List<String> extra = aes1.getAssemblyCategories().stream()
					.filter(x -> !pcCategories.contains(x))
					.collect(Collectors.toList());
			if (extra.size() > 0) {
				String msg = String.format(
						"Fatal error: Found %s in GRIDSS assembly files %s but not in input." +
								" GRIDSS does not support including assembly of reads included in INPUT files." +
								" If labels are used, they must match for both assembly and variant calling.",
						extra.get(0),
						aes1.getFile());
				throw new IllegalArgumentException(msg);
			}
			// Assembly duplicate
			for (int j = i + 1; j < aesList.size(); j++) {
				AssemblyEvidenceSource aes2 = aesList.get(j);
				List<String> common = aes1.getAssemblyCategories().stream()
						.filter(x -> aes2.getAssemblyCategories().contains(x))
						.collect(Collectors.toList());
				if (common.size() > 0) {
					String msg = String.format(
							"Fatal error: Found %s in GRIDSS assembly files %s and %s. " +
									"Each input file/label must be assembled in exactly one assembly.",
							common.get(0),
							aes1.getFile(),
							aes2.getFile());
					throw new IllegalArgumentException(msg);
				}
			}
		}
		// Missing in assembly
		List<String> allAssembledCategories = aesList.stream()
				.flatMap(aes -> aes.getAssemblyCategories().stream())
				.collect(Collectors.toList());
		for (int i = 0; i < pc.getCategoryCount(); i++) {
			String label = pc.getCategoryLabel(i);
			if (!allAssembledCategories.contains(label)) {
				String msg = String.format(
						"Fatal error: Missing assembly for %s. All input files must have a corresponding assembly.",
						label);
				throw new IllegalArgumentException(msg);
			}
		}
	}

	public List<String> getAssemblyCategories() {
		if (assembledCategories == null) {
			assembledCategories = getAssemblyCategories(getHeader());
		}
		return assembledCategories;
	}

	public static List<String> getAssemblyCategories(SAMFileHeader header) {
		List<String> categories = header.getComments().stream()
				.map(s -> s.replaceFirst("@CO	", ""))
				.filter(s -> s.startsWith(INPUT_CATEGORY_SAM_HEADER_PREFIX))
				.map(s -> s.replaceFirst(INPUT_CATEGORY_SAM_HEADER_PREFIX, ""))
				.collect(Collectors.toList());
		return categories;
	}
	/**
	 * Transforms assembly categories into their processing context offsets.
	 * @return category for each assembly support ordinal
	 */
	public int[] getAssemblyCategoryToProcessingContextCategoryLookup() {
		if (assemblyOrdinalToProcessingCategoryLookup == null) {
			List<String> ac = getAssemblyCategories();
			assemblyOrdinalToProcessingCategoryLookup = new int[ac.size()];
			for (int j = 0; j < ac.size(); j++) {
				String assemblyLabel = ac.get(j);
				for (int i = 0; i < getContext().getCategoryCount(); i++) {
					String pcLabel = getContext().getCategoryLabel(i);
					if (pcLabel.equals(assemblyLabel)) {
						assemblyOrdinalToProcessingCategoryLookup[j] = i;
						break;
					}
				}
			}
		}
		return assemblyOrdinalToProcessingCategoryLookup;
	}

	public SAMFileHeader getHeader() {
		if (this.header == null) {
			if (!getFile().exists()) {
				this.header = getContext().getBasicSamHeader();
				if (!this.header.getComments().stream().anyMatch(s -> s.contains(INPUT_CATEGORY_SAM_HEADER_PREFIX))) {
					for (int i = 0; i < getContext().getCategoryCount(); i++) {
						String category = getContext().getCategoryLabel(i);
						this.header.addComment(INPUT_CATEGORY_SAM_HEADER_PREFIX + category);
					}
				}
			} else {
				try (SamReader reader = getReader()) {
					// getReader() populates the header for us
				} catch (IOException e) {
					throw new RuntimeIOException(e);
				}
			}
		}
		return this.header;
	}
}
