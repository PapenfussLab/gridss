package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.configuration.AssemblyConfiguration;
import au.edu.wehi.idsv.configuration.VisualisationConfiguration;
import au.edu.wehi.idsv.debruijn.ReadErrorCorrector;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.SamTags;
import au.edu.wehi.idsv.util.FileHelper;
import au.edu.wehi.idsv.util.FilenameUtil;
import au.edu.wehi.idsv.util.IntervalUtil;
import au.edu.wehi.idsv.visualisation.AssemblyTelemetry.AssemblyChunkTelemetry;
import au.edu.wehi.idsv.visualisation.PositionalDeBruijnGraphTracker;
import com.google.common.collect.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Log;
import org.apache.commons.io.FileUtils;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

/**
 * Assemblies non-reference breakend contigs
 * 
 * @author Daniel Cameron
 *
 */
public class PositionalAssembler implements Iterator<SAMRecord> {
	private static final Log log = Log.getInstance(PositionalAssembler.class);
	private final ProcessingContext context;
	private final AssemblyEvidenceSource source;
	private final AssemblyIdGenerator assemblyNameGenerator;
	private final PeekingIterator<DirectedEvidence> it;
	private final BreakendDirection direction;
	private NonReferenceContigAssembler currentAssembler = null;
	private String currentContig = "";
	private AssemblyChunkTelemetry telemetry = null;
	private final IntervalBed excludedRegions;
	private final IntervalBed safetyRegions;
	private EvidenceTracker evidenceTracker = null;
	private boolean contigGeneratedSinceException = false;
	public PositionalAssembler(ProcessingContext context, AssemblyEvidenceSource source, AssemblyIdGenerator assemblyNameGenerator, Iterator<DirectedEvidence> backingIterator, BreakendDirection direction, IntervalBed excludedRegions, IntervalBed safetyRegions) {
		this.context = context;
		this.source = source;
		this.assemblyNameGenerator = assemblyNameGenerator;
		this.direction = direction;
		this.excludedRegions = excludedRegions;
		this.safetyRegions = safetyRegions;
		if (direction != null) {
			backingIterator = Iterators.filter(backingIterator, x -> x.getBreakendSummary() != null && x.getBreakendSummary().direction == this.direction);
		}
		this.it = Iterators.peekingIterator(backingIterator);
	}
	public PositionalAssembler(ProcessingContext context, AssemblyEvidenceSource source, AssemblyIdGenerator assemblyNameGenerator, Iterator<DirectedEvidence> it, IntervalBed excludedRegions, IntervalBed safetyRegions) {
		this(context, source, assemblyNameGenerator, it, null, excludedRegions, safetyRegions);
	}
	@Override
	public boolean hasNext() {
		ensureAssembler(Defaults.ATTEMPT_ASSEMBLY_RECOVERY, null);
		return currentAssembler != null && currentAssembler.hasNext();
	}
	@Override
	public SAMRecord next() {
		ensureAssembler(Defaults.ATTEMPT_ASSEMBLY_RECOVERY, null);
		SAMRecord r = currentAssembler.next();
		if (direction != null) {
			// force assembly direction to match the direction supplied
			r.setAttribute(SamTags.ASSEMBLY_DIRECTION, direction.toChar());
		}
		contigGeneratedSinceException = true;
		return r;
	}
	private void flushIfRequired() {
		if (currentAssembler != null && !currentAssembler.hasNext()) {
			closeCurrentAssembler();
		}
	}
	private void closeCurrentAssembler() {
		if (currentAssembler != null && currentAssembler.getExportTracker() != null) {
			try {
				currentAssembler.getExportTracker().close();
			} catch (IOException e) {
				log.debug(e);
			}
		}
		currentAssembler = null;
	}
	private Set<DirectedEvidence> getEvidenceInCurrentAssembler() {
		Set<DirectedEvidence> reloadRecoverySet = new HashSet<>();
		if (currentAssembler != null) {
			Set<KmerEvidence> ske = evidenceTracker.getTrackedEvidence();
			ske.addAll(currentAssembler.getEvidenceUntrackedButNotYetReturnedByIterator());
			reloadRecoverySet = ske.stream()
					.map(ke -> ke.evidence())
					.collect(Collectors.toSet());
		}
		return reloadRecoverySet;
	}
	private void ensureAssembler(boolean attemptRecovery, Set<DirectedEvidence> preload) {
		try {
			ensureAssembler(preload);
		} catch (AssemblyThresholdReachedException atre) {
			Set<DirectedEvidence> reloadRecoverySet = getEvidenceInCurrentAssembler();
			// really aggressive error correction of every read might simplify our assembly graph
			ReadErrorCorrector.errorCorrect(context.getAssemblyParameters().errorCorrection.k, context.getAssemblyParameters().errorCorrection.kmerErrorCorrectionMultiple / 2, context.getAssemblyParameters().errorCorrection.maxCorrectionsInKmer, context.getAssemblyParameters().errorCorrection.deduplicateReadKmers, reloadRecoverySet);
			Set<DirectedEvidence> downsampledRecoverySet = downsampleEvidenceInRegion(reloadRecoverySet, atre.getRange());
			// restart assembly using the downsampled set of reads
			closeCurrentAssembler();
			ensureAssembler(true, downsampledRecoverySet);
		} catch (AssertionError|Exception e) {
			if (contigGeneratedSinceException) {
				contigGeneratedSinceException = false;
				Set<DirectedEvidence> reloadRecoverySet = getEvidenceInCurrentAssembler();
				String msg = String.format("Error during assembly of chromosome %s (%d reads in graph). Attempting recovery by rebuilding assembly graph.", currentContig, reloadRecoverySet.size());
				log.info(e, msg);
				StringWriter sw = new StringWriter();
				e.printStackTrace(new PrintWriter(sw));
				log.info(sw);
				closeCurrentAssembler();
				ensureAssembler(attemptRecovery, reloadRecoverySet);
			} else {
				File packagedFile = null;
				try {
					packagedFile = packageMinimalAssemblyErrorReproductionData(getEvidenceInCurrentAssembler());
				} catch (Exception ex) {
					log.error(ex, "Error packaging assembly error reproduction data.");
				}
				String msg = "Error assembling " + currentContig + ". Please raise an issue at https://github.com/PapenfussLab/gridss/issues";
				if (packagedFile != null && packagedFile.exists()) {
					msg = msg + ". If your data can be shared publicly, please also include a minimal data set for reproducing the problem by attaching " + packagedFile.toString() + ".";
				}
				boolean unableToRecover = false;
				if (!attemptRecovery) {
					log.error(e, msg);
					throw e;
				} else {
					try {
						if (it.hasNext()) {
							msg = String.format("%s. Attempting recovery by resuming assembly at %s:%d",
									msg,
									context.getReference().getSequenceDictionary().getSequence(it.peek().getBreakendSummary().referenceIndex).getSequenceName(),
									it.peek().getBreakendSummary().start);
						}
					} catch (AssertionError | Exception nested) {
						log.error(nested, "Assembly recovery attempt failed due to exception thrown by underlying iterator");
						unableToRecover = true;
					}
					log.error(e, msg);
				}
				if (unableToRecover) {
					throw e;
				} else {
					closeCurrentAssembler();
					// discard the reads that we already in the graph have and continue
					// check config as to whether we want to recover from back-to-back failures
					ensureAssembler(context.getConfig().getAssembly().recoverAfterError, null);
				}
			}
		}
	}

	private Set<DirectedEvidence> downsampleEvidenceInRegion(Set<DirectedEvidence> preload, Range<Integer> range) {
		float scRate = context.getAssemblyParameters().downsampling.densityDownsampleRateClippedReads;
		float rpRate = context.getAssemblyParameters().downsampling.densityDownsampleRateDiscordantReads;
		Random rng = new Random(range.lowerEndpoint());
		Set<DirectedEvidence> downsampledPreload = preload.stream().filter(de -> {
			float downsampleThreshold = scRate;
			boolean overlapsRanges = IntervalUtil.overlapsClosed(
					range.lowerEndpoint(), range.upperEndpoint(),
					de.getUnderlyingSAMRecord().getUnclippedStart(), de.getUnderlyingSAMRecord().getUnclippedEnd());
			if (de instanceof NonReferenceReadPair) {
				downsampleThreshold = rpRate;
				NonReferenceReadPair nrrp = (NonReferenceReadPair) de;
				KmerEvidence ke = KmerEvidence.create(context.getAssemblyParameters().k, nrrp);
				overlapsRanges |= IntervalUtil.overlapsClosed(
						range.lowerEndpoint(), range.upperEndpoint(),
						ke.startPosition(), ke.endPosition());
			}
			return !overlapsRanges | rng.nextFloat() < downsampleThreshold;
		}).collect(Collectors.toSet());
		if (downsampledPreload.size() == preload.size()) {
			log.warn(String.format("Forcing downsampling around %s:%d-%d since probabilistic downsampling failed.", currentContig, range.lowerEndpoint(), range.upperEndpoint()));
			// force an overall downsampling if we didn't actually remove anything
			downsampledPreload = downsampledPreload.stream()
					.skip(1)
					.filter(de -> rng.nextFloat() < 0.95)
					.collect(Collectors.toSet());
		}
		return downsampledPreload;
	}


	private static AtomicInteger errorPackagesCreated = new AtomicInteger(0);
	private File packageMinimalAssemblyErrorReproductionData(Set<DirectedEvidence> evidenceInCurrentAssembler) throws IOException {
		FileSystemContext fsc = context.getFileSystemContext();
		File asswd = fsc.getWorkingDirectory(source.getFile());
		int fatalErrorNumber = errorPackagesCreated.incrementAndGet();
		File outFile = new File(asswd, String.format("gridss_minimal_reproduction_data_for_error_%d.zip", fatalErrorNumber));
		int maxPackages = context.getConfig().getAssembly().maximumReproductionExportPackages;
		if (fatalErrorNumber > maxPackages) {
			log.warn("Already exported " + maxPackages + " minimal reproduction data exports. Not generating " + outFile);
			return null;
		}
		File workingLocation = new File(asswd, outFile.getName() + ".error_package");
		if (!asswd.exists()) {
			asswd.mkdir();
		}
		workingLocation.mkdir();
		// Copy config file
		if (context.getConfig().getSourceConfigurationFile() != null) {
			FileUtils.copyFileToDirectory(context.getConfig().getSourceConfigurationFile(), workingLocation);
		}
		Set<SAMEvidenceSource> sources = evidenceInCurrentAssembler.stream()
				.map(de -> (SAMEvidenceSource)de.getEvidenceSource())
				.collect(Collectors.toSet());
		for (SAMEvidenceSource source : sources) {
			File samFile = source.getFile();
			File samInDir = fsc.getIntermediateDirectory(samFile);
			File samOutDir = new File(workingLocation, samInDir.getName());
			samOutDir.mkdir();
			// Copy metrics
			FileUtils.copyFileToDirectory(fsc.getCigarMetrics(samFile), samOutDir);
			FileUtils.copyFileToDirectory(fsc.getCoverageBlacklistBed(samFile), samOutDir);
			FileUtils.copyFileToDirectory(fsc.getIdsvMetrics(samFile), samOutDir);
			FileUtils.copyFileToDirectory(fsc.getInsertSizeMetrics(samFile), samOutDir);
			FileUtils.copyFileToDirectory(fsc.getMapqMetrics(samFile), samOutDir);
			FileUtils.copyFileToDirectory(fsc.getTagMetrics(samFile), samOutDir);

			// Copy minimal set of reads
			File outSam = new File(samOutDir, fsc.getSVBam(samFile).getName());
			SAMFileHeader header;
			try (SamReader reader = context.getSamReader(source.getFile())) {
				header = reader.getFileHeader();
			}
			// TODO: preserve original BAM ordering?
			try (SAMFileWriter writer = context.getSamFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, outSam)) {
				evidenceInCurrentAssembler.stream()
						.filter(de -> de.getEvidenceSource() == source)
						.forEach(a -> writer.addAlignment(a.getUnderlyingSAMRecord()));
			}
		}
		// reference genome - regions not overlapping support are N masked to improve compression
		final int paddingMargin = Math.min(2000, 2 * source.getMaxConcordantFragmentSize());
		ReferenceLookup rl = context.getReference();
		File ref = new File(workingLocation, "masked_ref.fa");
		try (FileOutputStream fso = new FileOutputStream(ref)) {
			for (SAMSequenceRecord seq : rl.getSequenceDictionary().getSequences()) {
				final int refIndex = seq.getSequenceIndex();
				fso.write(new byte[] { '>'});
				fso.write(seq.getSequenceName().getBytes(StandardCharsets.UTF_8));
				fso.write(new byte[] { '\n'});
				RangeSet<Integer> rs = TreeRangeSet.create(Streams.concat(
					// local side of evidence
					evidenceInCurrentAssembler.stream()
						.map(de -> de.getBreakendSummary())
						.filter(bs -> bs.referenceIndex == refIndex)
						.map(bs -> Range.closed(bs.start - paddingMargin, bs.end + paddingMargin)),
					// remote side of evidence
					evidenceInCurrentAssembler.stream()
						.map(de -> de.getBreakendSummary())
						.filter(bs -> bs instanceof BreakpointSummary)
						.map(bs -> (BreakpointSummary)bs)
						.filter(bs -> bs.referenceIndex2 == refIndex)
						.map(bs -> Range.closed(bs.start2 - paddingMargin, bs.end2 + paddingMargin)))
					.collect(Collectors.toList()));
				byte[] bases = rl.getSequence(seq.getSequenceName()).getBases().clone();
				if (rs.isEmpty()) {
					Arrays.fill(bases, (byte)'N');
				} else {
					for (Range<Integer> nRange : rs.complement().asRanges()) {
						int lowerBound = Math.max(0, nRange.hasLowerBound() ? nRange.lowerEndpoint() : 0);
						int upperBound = Math.min(bases.length, nRange.hasUpperBound() ? nRange.upperEndpoint() : bases.length);
						upperBound = Math.max(lowerBound, upperBound);
						Arrays.fill(bases, lowerBound, upperBound, (byte) 'N');
					}
				}
				fso.write(bases);
				fso.write(new byte[] { '\n'});
			}
		}
		FileHelper.zipDirectory(outFile, workingLocation);
		FileUtils.deleteDirectory(workingLocation);
		return outFile;
	}

	private void ensureAssembler(Set<DirectedEvidence> preload) {
		flushIfRequired();
		while ((currentAssembler == null || !currentAssembler.hasNext()) && ((preload != null && !preload.isEmpty()) || it.hasNext())) {
			// traverse contigs until we find one that has an assembly to call
			currentAssembler = createAssembler(preload);
			preload = null;
			flushIfRequired();
		}
	}
	private NonReferenceContigAssembler createAssembler(Set<DirectedEvidence> preload) {
		AssemblyConfiguration ap = context.getAssemblyParameters();
		int maxKmerSupportIntervalWidth = source.getMaxConcordantFragmentSize() - source.getMinConcordantFragmentSize() + 1; 		
		int maxReadLength = source.getMaxReadLength();
		int k = ap.k;
		int maxEvidenceSupportIntervalWidth = maxKmerSupportIntervalWidth + maxReadLength - k + 2;
		int maxPathLength = ap.positional.maxPathLengthInBases(maxReadLength);
		int anchorAssemblyLength = ap.anchorLength;
		PeekingIterator<DirectedEvidence> inputIterator = it;
		if (preload != null && preload.size() > 0) {
			ArrayList<DirectedEvidence> list = Lists.newArrayList(preload);
			list.sort(DirectedEvidenceOrder.ByNatural);
			inputIterator = Iterators.peekingIterator(Iterators.concat(list.iterator(), inputIterator));
		}
		int referenceIndex = inputIterator.peek().getBreakendSummary().referenceIndex;
		int firstPosition = inputIterator.peek().getBreakendSummary().start;
		currentContig = context.getDictionary().getSequence(referenceIndex).getSequenceName();
		ReferenceIndexIterator evidenceIt = new ReferenceIndexIterator(inputIterator, referenceIndex);
		evidenceTracker = new EvidenceTracker();
		SupportNodeIterator supportIt = new SupportNodeIterator(k, evidenceIt, Math.max(2 * source.getMaxReadLength(), source.getMaxConcordantFragmentSize()), evidenceTracker, ap.includePairAnchors, ap.pairAnchorMismatchIgnoreEndBases);
		AggregateNodeIterator agIt = new AggregateNodeIterator(supportIt);
		Iterator<KmerNode> knIt = agIt;
		if (Defaults.SANITY_CHECK_ASSEMBLY_GRAPH) {
			knIt = evidenceTracker.new AggregateNodeAssertionInterceptor(knIt);
		}
		PathNodeIterator pathNodeIt = new PathNodeIterator(knIt, maxPathLength, k); 
		Iterator<KmerPathNode> pnIt = pathNodeIt;
		if (Defaults.SANITY_CHECK_ASSEMBLY_GRAPH) {
			pnIt = evidenceTracker.new PathNodeAssertionInterceptor(pnIt, "PathNodeIterator");
		}
		currentAssembler = new NonReferenceContigAssembler(pnIt, referenceIndex, maxEvidenceSupportIntervalWidth, anchorAssemblyLength, k, source, assemblyNameGenerator, evidenceTracker, currentContig, BreakendDirection.Forward, excludedRegions, safetyRegions);
		VisualisationConfiguration vis = context.getConfig().getVisualisation();
		if (vis.assemblyProgress) {
			String filename = String.format("positional-%s_%d-%s.csv", context.getDictionary().getSequence(referenceIndex).getSequenceName(), firstPosition, direction);
			filename = FilenameUtil.stripInvalidFilenameCharacters(filename);
			File file = new File(vis.directory, filename);
			PositionalDeBruijnGraphTracker exportTracker;
			try {
				exportTracker = new PositionalDeBruijnGraphTracker(file, supportIt, agIt, pathNodeIt, null, evidenceTracker, currentAssembler);
				exportTracker.writeHeader();
				currentAssembler.setExportTracker(exportTracker);
			} catch (IOException e) {
				log.debug(e);
			}
		}
		currentAssembler.setTelemetry(getTelemetry());
		return currentAssembler;
	}
	public AssemblyChunkTelemetry getTelemetry() {
		return telemetry;
	}
	public void setTelemetry(AssemblyChunkTelemetry assemblyChunkTelemetry) {
		this.telemetry = assemblyChunkTelemetry;
	}
	private static class ReferenceIndexIterator implements PeekingIterator<DirectedEvidence> {
		private final PeekingIterator<DirectedEvidence> it;
		private final int referenceIndex;
		public ReferenceIndexIterator(PeekingIterator<DirectedEvidence> it, int referenceIndex) {
			this.it = it;
			this.referenceIndex = referenceIndex;
		}
		@Override
		public boolean hasNext() {
			return it.hasNext() && it.peek().getBreakendSummary().referenceIndex == referenceIndex;
		}

		@Override
		public DirectedEvidence next() {
			if (!hasNext()) throw new NoSuchElementException();
			return it.next();
		}

		@Override
		public DirectedEvidence peek() {
			if (hasNext()) return it.peek();
			throw new NoSuchElementException();
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}

	}
}
