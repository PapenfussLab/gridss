package au.edu.wehi.idsv;

import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.configuration.GridssConfiguration;
import au.edu.wehi.idsv.configuration.SoftClipConfiguration;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.sam.SAMFileUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.*;
import au.edu.wehi.idsv.validation.OrderAssertingIterator;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import gridss.ComputeSamTags;
import gridss.ExtractSVReads;
import gridss.SoftClipsToSplitReads;
import gridss.analysis.CollectGridssMetrics;
import gridss.cmdline.CommandLineProgramHelper;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.Path;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Structural variation evidence based on read pairs from a single SAM/BAM.  
 * @author Daniel Cameron
 *
 */
public class SAMEvidenceSource extends EvidenceSource {
	private static final Log log = Log.getInstance(SAMEvidenceSource.class);
	private final int sourceCategory;
	private final Integer rpcMinFragmentSize;
	private final Integer rpcMaxFragmentSize;
	private final Double rpcConcordantPercentage;
	private IdsvSamFileMetrics metrics;
	private ReadPairConcordanceCalculator rpcc;
	public enum EvidenceSortOrder {
		SAMRecordStartPosition,
		EvidenceStartPosition
	}
	public SAMEvidenceSource(ProcessingContext processContext, File file, File nameSorted, int sourceCategory) {
		this(processContext, file, nameSorted, sourceCategory, null, null, null);
	}
	public SAMEvidenceSource(ProcessingContext processContext, File file, File nameSorted, int sourceCategory, int minFragmentSize, int maxFragmentSize) {
		this(processContext, file, nameSorted, sourceCategory, minFragmentSize, maxFragmentSize, null);
	}
	public SAMEvidenceSource(ProcessingContext processContext, File file, File nameSorted, int sourceCategory, double concordantPercentage) {
		this(processContext, file, nameSorted, sourceCategory, null, null, concordantPercentage);
	}
	protected SAMEvidenceSource(ProcessingContext processContext, File file, File nameSorted,
			int sourceCategory, Integer rpcMinFragmentSize, Integer rpcMaxFragmentSize, Double rpcConcordantPercentage) {
		super(processContext, file, nameSorted);
		this.sourceCategory = sourceCategory;
		this.rpcMinFragmentSize = rpcMinFragmentSize;
		this.rpcMaxFragmentSize = rpcMaxFragmentSize;
		this.rpcConcordantPercentage = rpcConcordantPercentage;
	}
	public IdsvSamFileMetrics getMetrics() {
		if (metrics == null) {
			ensureMetrics();
		}
		return metrics;
	}
	
	/**
	 * Grouping category of source. For tumour/normal analysis, normal is considered group 0, tumour group 1
	 * @return
	 */
	public int getSourceCategory() {
		return sourceCategory;
	}
	public synchronized void ensureMetrics() {
		if (metrics == null) {
			File idsvFile = getContext().getFileSystemContext().getIdsvMetrics(getFile());
			File cigarFile = getContext().getFileSystemContext().getCigarMetrics(getFile());
			File mapqFile = getContext().getFileSystemContext().getMapqMetrics(getFile());
			if (!idsvFile.exists() || !cigarFile.exists() || !mapqFile.exists()) {
				log.info("Calculating metrics for " + getFile().getAbsolutePath());
				CommandLineProgramHelper cmd = new CommandLineProgramHelper(new CollectGridssMetrics());
				cmd.addArg("INPUT", getFile().getPath());
				cmd.addArg("OUTPUT", getContext().getFileSystemContext().getMetricsPrefix(getFile()).getPath());
				cmd.addArg("THRESHOLD_COVERAGE", getContext().getConfig().maxCoverage);
				cmd.addArg("FILE_EXTENSION", "null");
				cmd.addArg("GRIDSS_PROGRAM", "CollectCigarMetrics");
				cmd.addArg("GRIDSS_PROGRAM", "CollectMapqMetrics");
				cmd.addArg("GRIDSS_PROGRAM", "CollectTagMetrics");
				cmd.addArg("GRIDSS_PROGRAM", "CollectIdsvMetrics");
				cmd.addArg("GRIDSS_PROGRAM", "ReportThresholdCoverage");
				cmd.addArg("PROGRAM", "null");
				if (!knownSingleEnded()) {
					// Don't run CollectInsertSizeMetrics
					cmd.addArg("PROGRAM","CollectInsertSizeMetrics");
				} else {
					// The CollectMultipleMetrics super class complains if no PROGRAM set so
					// we'll just collect some stuff that is useful, but we don't actually
					// use yet
					cmd.addArg("PROGRAM","CollectAlignmentSummaryMetrics");
				}
				if (getContext().getCalculateMetricsRecordCount() < Integer.MAX_VALUE) {
					cmd.addArg("STOP_AFTER", getContext().getCalculateMetricsRecordCount());
				}
				execute(cmd);
			}
			metrics = new IdsvSamFileMetrics(getContext(), getFile(), knownSingleEnded());
		}
	}
	protected void execute(CommandLineProgramHelper cmd) {
		if (cmd.getProgram() instanceof ReferenceCommandLineProgram) {
			((ReferenceCommandLineProgram) cmd.getProgram()).setReference(getContext().getReference());
		}
		if (getContext().getCommandLineProgram() == null) {
			cmd.addArg("REFERENCE_SEQUENCE", getContext().getReferenceFile());
			cmd.addArg("TMP_DIR", getContext().getFileSystemContext().getTemporaryDirectory());
			cmd.addArg("MAX_RECORDS_IN_RAM", getContext().getFileSystemContext().getMaxBufferedRecordsPerFile());
		} else {
			cmd.setCommonArgs(getContext().getCommandLineProgram());
		}
		int result = cmd.run();
		if (result != 0) {
			String msg = "Unable to execute " + cmd.getClass().getName() + " for " + getFile();
			log.error(msg);
			if (getContext().getConfig().terminateOnFirstError) {
				System.exit(result);
			}
			throw new RuntimeException(msg);
		}
	}
	public synchronized void ensureExtracted() throws IOException {
		File svFile = getContext().getFileSystemContext().getSVBam(getFile());
		File extractedFile = FileSystemContext.getWorkingFileFor(svFile, "gridss.tmp.extracted.");
		File querysortedFile = FileSystemContext.getWorkingFileFor(svFile, "gridss.tmp.querysorted.");
		File taggedFile = FileSystemContext.getWorkingFileFor(svFile, "gridss.tmp.tagged.");
		File withsplitreadsFile = FileSystemContext.getWorkingFileFor(svFile, "gridss.tmp.splitreads.");
		ensureMetrics();
		// Regenerate from from the intermediate file furtherest through the pipeline
		// extract -> query sort -> tag -> split read -> back to coordinate sorted
		// We want to tag before generating split reads so all splits are guaranteed to
		// have the same tags
		if (!svFile.exists()) {
			if (!withsplitreadsFile.exists()) {
				if (!taggedFile.exists()) {
					if (!querysortedFile.exists()) {
						if (!extractedFile.exists()) {
							log.info("Extracting SV reads from " + getFile().getAbsolutePath());
							File in = getFile(SortOrder.queryname);
							if (in == null || !in.exists()) {
								in = getFile();
							}
							CommandLineProgramHelper cmd = new CommandLineProgramHelper(new ExtractSVReads());
							cmd.addArg("INPUT", in.getPath());
							cmd.addArg("OUTPUT", extractedFile.getPath());
							cmd.addArg("UNMAPPED_READS", "false"); // saves intermediate file space
							cmd.addArg("MIN_CLIP_LENGTH", getContext().getConfig().getSoftClip().minLength);
							cmd.addArg("INSERT_SIZE_METRICS", getContext().getFileSystemContext().getInsertSizeMetrics(getFile()));
							// Picard tools does not mark duplicates correctly. We need to keep them so we can
							// fix the duplicate marking in ComputeSamTags
							cmd.addArg("INCLUDE_DUPLICATES", "true");
							if (rpcMinFragmentSize != null) cmd.addArg("READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE", rpcMinFragmentSize);
							if (rpcMaxFragmentSize != null) cmd.addArg("READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE", rpcMaxFragmentSize);
							if (rpcConcordantPercentage != null) cmd.addArg("READ_PAIR_CONCORDANT_PERCENT", rpcConcordantPercentage);
							execute(cmd);
						}
						SAMFileUtil.sort(getContext().getFileSystemContext(), extractedFile, querysortedFile, SortOrder.queryname);
						if (gridss.Defaults.DELETE_TEMPORARY_FILES) {
							FileHelper.delete(extractedFile, true);
						}
					}
					log.info("Computing SAM tags for " + svFile);
					CommandLineProgramHelper cmd = new CommandLineProgramHelper(new ComputeSamTags());
					cmd.addArg("INPUT", querysortedFile.getPath());
					cmd.addArg("OUTPUT", taggedFile.getPath());
					execute(cmd);
					if (gridss.Defaults.DELETE_TEMPORARY_FILES) {
						FileHelper.delete(querysortedFile, true);
					}
				}
				log.info("Identifying split reads for " + getFile().getAbsolutePath());
				SoftClipsToSplitReads program = new SoftClipsToSplitReads();
				CommandLineProgramHelper cmd = new CommandLineProgramHelper(program);
				cmd.addArg("WORKER_THREADS", getProcessContext().getWorkerThreadCount());
				cmd.addArg("INPUT", taggedFile.getPath());
				cmd.addArg("OUTPUT", withsplitreadsFile.getPath());
				cmd.addArg("REALIGN_EXISTING_SPLIT_READS", Boolean.toString(getContext().getConfig().getSoftClip().realignSplitReads));
						// realignment.* not soft-clip
						//"MIN_CLIP_LENGTH=" + getContext().getConfig().
						//"MIN_CLIP_QUAL=" + getContext().getConfig().getSoftClip().minAverageQual);
				program.setReference(getProcessContext().getReference());
				program.setFileSystemContext(getProcessContext().getFileSystemContext());
				execute(cmd);
				if (gridss.Defaults.DELETE_TEMPORARY_FILES) {
					FileHelper.delete(taggedFile, true);
				}
			}
			SAMFileUtil.sort(getContext().getFileSystemContext(), withsplitreadsFile, svFile, SortOrder.coordinate);
			if (gridss.Defaults.DELETE_TEMPORARY_FILES) {
				FileHelper.delete(withsplitreadsFile, true);
			}
		}
		if (gridss.Defaults.DELETE_TEMPORARY_FILES) {
			FileHelper.delete(extractedFile, true);
			FileHelper.delete(querysortedFile, true);
			FileHelper.delete(taggedFile, true);
			FileHelper.delete(withsplitreadsFile, true);
		}
	}
	public CloseableIterator<DirectedEvidence> iterator(final QueryInterval[] intervals, EvidenceSortOrder eso) {
		SamReader reader = getReader();
		// expand query bounds as the alignment for a discordant read pair could fall before or after the breakend interval we are extracting
		QueryInterval[] expandedIntervals = QueryIntervalUtil.padIntervals(getContext().getDictionary(), intervals, getMaxConcordantFragmentSize() + 1);
		// ignore blacklisted regions
		IntervalBed queryInterval = new IntervalBed(getContext().getLinear(), expandedIntervals);
		queryInterval.remove(getBlacklistedRegions());
		CloseableIterator it = tryOpenReader(reader, queryInterval.asQueryInterval());
		if (Defaults.SANITY_CHECK_DUMP_ITERATORS) {
			it = new AutoClosingIterator<>(new DebugSpammingIterator<>(it, "SAMEvidenceSource.iterator().rawiterator"));
		}
		Iterator<DirectedEvidence> eit = asEvidence(it, eso);
		eit = Iterators.filter(eit, e -> QueryIntervalUtil.overlaps(intervals, e.getBreakendSummary()));
		if (Defaults.SANITY_CHECK_DUMP_ITERATORS) {
			eit = new AutoClosingIterator<>(new DebugSpammingIterator<>(eit, "SAMEvidenceSource.iterator(QueryInterval[]).filter"));
		}
		return new AutoClosingIterator<>(eit, reader, it);
	}
	/**
	 * Attempts to open a new iterator.
	 * 
	 * As htsjdk performing memory mapping internally, many open/close cycles can results in an exhaustion
	 * of file handles. 
	 * 
	 * @param reader
	 * @param intervals
	 * @return
	 */
	private SAMRecordIterator tryOpenReader(SamReader reader, QueryInterval[] intervals) {
		SAMRecordIterator it = null;
		try {
			it = reader.queryOverlapping(intervals);
		} catch (Exception e) {
			log.debug("Attempting to recover from query failure: ", e);
			System.gc();
			System.runFinalization();
			System.gc();
			it = reader.queryOverlapping(intervals);
			log.debug("Recovery successful");
		}
		return it;
	}
	public CloseableIterator<DirectedEvidence> iterator(EvidenceSortOrder eso) {
		SamReader reader = getReader();
		SAMRecordIterator it = reader.iterator();
		it.assertSorted(SortOrder.coordinate);
		Iterator<DirectedEvidence> eit = asEvidence(it, eso);
		if (Defaults.SANITY_CHECK_DUMP_ITERATORS) {
			eit = new AutoClosingIterator<>(new DebugSpammingIterator<>(eit, "SAMEvidenceSource.iterator()"));
		}
		return new AutoClosingIterator<>(eit, reader, it);
	}
	protected SamReader getReader() {
		File svFile = getSVFile();
		SamReader reader = getProcessContext().getSamReader(svFile.exists() ? svFile : getFile());
		return reader;
	}

	public File getSVFile() {
		if (getFile() == null) {
			return null;
		}
		return getContext().getFileSystemContext().getSVBam(getFile());
	}
	public void assertPreprocessingComplete() {
		File svFile = getSVFile();
		if (svFile != null && !svFile.exists()) {
			throw new IllegalStateException(String.format("Missing required file %s. See GRIDSS pipeline examples and documentation.", svFile));
		}
	}
	private Iterator<DirectedEvidence> asEvidence(Iterator<SAMRecord> it, EvidenceSortOrder eso) {
		it = new BufferedIterator<>(it, 2); // TODO: remove when https://github.com/samtools/htsjdk/issues/760 is resolved
		it = Iterators.filter(it, r -> !shouldFilterPreTransform(r));
		if (Defaults.SANITY_CHECK_DUMP_ITERATORS) {
			it = new AutoClosingIterator<>(new DebugSpammingIterator<>(it, "SAMEvidenceSource.shouldFilterPreTransform"));
		}
		it = Iterators.transform(it, r -> transform(r));
		if (Defaults.SANITY_CHECK_DUMP_ITERATORS) {
			it = new AutoClosingIterator<>(new DebugSpammingIterator<>(it, "SAMEvidenceSource.transform"));
		}
		it = Iterators.filter(it, r -> !shouldFilter(r));
		if (Defaults.SANITY_CHECK_DUMP_ITERATORS) {
			it = new AutoClosingIterator<>(new DebugSpammingIterator<>(it, "SAMEvidenceSource.shouldFilter(SAMRecord)"));
		}
		Iterator<DirectedEvidence> eit = new DirectedEvidenceIterator(it, this, minIndelSize());
		eit = Iterators.filter(eit, e -> !shouldFilter(e));
		if (Defaults.SANITY_CHECK_DUMP_ITERATORS) {
			it = new AutoClosingIterator<>(new DebugSpammingIterator<>(it, "SAMEvidenceSource.shouldFilter(DirectedEvidence)"));
		}
		switch (eso) {
			case SAMRecordStartPosition:
				// already sorted by coordinate
				break;
			case EvidenceStartPosition:
				eit = new DirectEvidenceWindowedSortingIterator<DirectedEvidence>(getContext(), getSortWindowSize(), eit);
				if (Defaults.SANITY_CHECK_ITERATORS) {
					// Can't enforce pairing as there may actually be duplicates depending on how multi-mapping alignment was performed
					//eit = new PairedEvidenceTracker<DirectedEvidence>(getFile().getName(), eit);
					eit = new OrderAssertingIterator<DirectedEvidence>(eit, DirectedEvidenceOrder.ByNatural);
				}
				break;
			default:
				throw new IllegalArgumentException("Sort order must be specified");
		}
		if (Defaults.SANITY_CHECK_DUMP_ITERATORS) {
			eit = new DebugSpammingIterator<>(eit, "asEvidence");
		}
		return eit;
	}
	private static float average(byte[] values) {
		float total = 0;
		for (byte b : values) {
			total += b;
		}
		return total / values.length;
	}
	public SAMRecord transform(SAMRecord r) {
		SAMRecordUtil.lowMapqToUnmapped(r, getContext().getConfig().minMapq);
		// Converts overlaps of blacklisted regions to unmapped
		if (!r.getReadUnmappedFlag()) {
			if (getBlacklistedRegions().overlaps(r.getReferenceIndex(), r.getAlignmentStart(), r.getAlignmentEnd())) {
				r.setReadUnmappedFlag(true);
			}
		}
		if (r.getReadPairedFlag() && !r.getMateUnmappedFlag()) {
			int mateRef = r.getMateReferenceIndex();
			int mateStart = r.getMateAlignmentStart();
			int mateEnd = mateStart;
			Cigar mateCigar = SAMRecordUtil.getCachedMateCigar(r);
			if (mateCigar != null) {
				mateEnd += mateCigar.getReferenceLength() - 1;
			}
			if (getBlacklistedRegions().overlaps(mateRef, mateStart, mateEnd)) {
				r.setMateUnmappedFlag(true);
			}
		}
		if (r.getAttribute(SAMTag.SA.name()) != null) {
			SAMSequenceDictionary dict = getContext().getDictionary();
			// Need to keep track of original SA tag as if we unmap the primary alignment
			// the supplementary alignment scoring will be inconsistent since it is based
			// on the length of the primary alignment soft clip.
			r.setTransientAttribute("OSA", r.getStringAttribute(SAMTag.SA.name()));
			r.setAttribute(SAMTag.SA.name(), ChimericAlignment.getChimericAlignments(r).stream()
					.filter(ca -> isInReference(r, ca, dict))
					.filter(ca -> !getBlacklistedRegions().overlaps(
							dict.getSequence(ca.rname).getSequenceIndex(),
							ca.pos,
							ca.pos + ca.cigar.getReferenceLength() - 1))
					.map(ca -> ca.toString())
					.collect(Collectors.joining(";")));
		}
		return r;
	}

	private static boolean isInReference(SAMRecord r, ChimericAlignment ca, SAMSequenceDictionary dict) {
		SAMSequenceRecord seq = dict.getSequence(ca.rname);
		if (seq == null) {
			if (!MessageThrottler.Current.shouldSupress(log, "SA contig")) {
				log.error(String.format("Read %s contains split read (SA tag) reference to \"%s\" which does not exist in the reference. Ignoring"), r.getReadName(), ca.rname);
			}
		}
		return seq != null;
	}

	/**
	 * Fast filtering logic to reduce the number of records that need to be converted
	 * from SAMRecords to DirectedEvidence.
	 * @param r
	 * @return
	 */
	public boolean shouldFilterPreTransform(SAMRecord r) {
		if (r == null) {
			return true;
		}
		if (r.getReadUnmappedFlag() || r.getMappingQuality() < getContext().getConfig().minMapq) {
			return true;
		}
		if (getContext().isFilterDuplicates() && r.getDuplicateReadFlag()) {
			return true;
		}
		if (!isIndelOrClipped(r)) {
			if (!r.getReadPairedFlag() || (rpcc != null && rpcc.isConcordant(r))) {
				return true;
			}
		}
		return false;
	}
	private boolean isIndelOrClipped(SAMRecord r) {
		for (CigarElement ce : r.getCigar()) {
			switch (ce.getOperator()) {
				case S:
				case H:
				case D:
				case I:
				case N:
					return true;
			}
		}
		return false;
	}
	public boolean shouldFilter(SAMRecord r) {
		if (r == null) {
			return true;
		}
		if (r.getReadUnmappedFlag()) {
			return true;
		}
		if (getContext().isFilterDuplicates() && r.getDuplicateReadFlag()) {
			return true;
		}
		if (CigarUtil.widthOfImprecision(r.getCigar()) == 0) {
			if (r.getAlignmentStart() < 1) {
				return true;
			}
			if (r.getAlignmentEnd() > getProcessContext().getReference().getSequenceDictionary().getSequence(r.getReferenceIndex()).getSequenceLength()) {
				return true;
			}
		}
		return false;
	}
	private IntervalBed blacklist = null;
	public IntervalBed getBlacklistedRegions() {
		if (blacklist == null) {
			File coverageBlacklist = getContext().getFileSystemContext().getCoverageBlacklistBed(getFile());
			if (!coverageBlacklist.exists()) {
				// Fall back to generic blacklist if we haven't calculated a coverage blacklist yet 
				blacklist = getContext().getBlacklistedRegions();
			} else {
				try {
					blacklist = IntervalBed.merge(getContext().getLinear(), ImmutableList.of(
							getContext().getBlacklistedRegions(),
							new IntervalBed(getContext().getLinear(), coverageBlacklist)
					));
				} catch (IOException e) {
					log.error(e);
					blacklist = getContext().getBlacklistedRegions();
				}
			}
			if (blacklist == null) {
				blacklist = new IntervalBed(getContext().getLinear());
			}
		}
		return blacklist;
	}
	// Exposed mostly for testing purposes
	protected void setBlacklistedRegions(IntervalBed blacklist) {
		this.blacklist = blacklist;
	}
	private int minIndelSize() {
		return Math.min(getContext().getConfig().getSoftClip().minLength, getContext().getVariantCallingParameters().minSize);
	}
	public boolean shouldFilter(DirectedEvidence e) {
		BreakendSummary bs = e.getBreakendSummary();
		if (getBlacklistedRegions().overlaps(bs.referenceIndex, bs.start - 1, bs.end + 1)) {
			return true;
		}
		GridssConfiguration config = getContext().getConfig();
		if (e instanceof SingleReadEvidence) {
			if (((SingleReadEvidence) e).isReference()) {
				return true;
			}
		}
		if (e instanceof SoftClipEvidence) {
			SoftClipEvidence sce = (SoftClipEvidence) e;
			SoftClipConfiguration scc = config.getSoftClip();
			if (sce.getBreakendSequence().length < scc.minLength) return true;
			if (sce.getBreakendQuality() != null && average(sce.getBreakendQuality()) < scc.minAverageQual) return true;
			if (config.adapters.isAdapterSoftClip(sce)) return true;
			if (!AssemblyAttributes.isAssembly(sce.getSAMRecord())) {
				// TODO: symmetrical identity and entropy filters on both sides
				if(sce.getEvidenceSource().getFile().getPath().toLowerCase().endsWith(".cram")){
					// in CRAM files, the MD and NM tags are not always present
					htsjdk.samtools.util.SequenceUtil.calculateMdAndNmTags(sce.getSAMRecord(),
							getContext().getReferenceSequence(sce.getSAMRecord().getContig()).getBases(), true, true);
				}
				if (SAMRecordUtil.getAlignedIdentity(sce.getSAMRecord()) < scc.minAnchorIdentity) return true;
				if (SAMRecordUtil.alignedEntropy(sce.getSAMRecord()) < config.minAnchorShannonEntropy) return true;
			}
		}
		if (e instanceof IndelEvidence) {
			// Currently filtering small indels in evidence iterator itself
		}
		if (e instanceof DirectedBreakpoint) {
			BreakpointSummary bp = (BreakpointSummary)e.getBreakendSummary();
			if (getBlacklistedRegions().overlaps(bp.referenceIndex2, bp.start2 - 1, bp.end2 + 1)) {
				return true;
			}
			// Still do assembly - leave the filtering to the variant calling
			//if (bp.getEventSize() != null && bp.getEventSize() < config.getVariantCalling().minSize) {
			//	return true;
			//}
		}
		return false;
	}
	protected int getSortWindowSize() {
		// Soft clip:
		// worst case: forward with small clip followed by large homologous clip 
		// MMMMMMMMMM>
		//          <<<<<<<<<<M with microhomology the full length of the read
		int softClipSortWindowSize = 2 * Math.max(getMaxReadMappedLength(), getMaxReadLength()) + 1;
		// Read pair:
		// worst case scenario is a fully mapped forward read followed by a large soft-clipped backward read at max distance
		// MDDDDDDDDDDDDDDDDDMSSSSSSSSSSSSSSS>
		// ^--mapped length--^
		//                  ^--------max concordant fragment size-----^
		//                  |-----       breakend call            ----|
		//                                                   <SSSSSSSSSM
		//                   ^--------max concordant fragment size-----^
		// ^ alignment start                                           ^ alignment start
		int readPairSortWindowSize = getMaxConcordantFragmentSize() + getMaxReadMappedLength() + getMaxReadLength() + 1;
		int windowSize = Math.max(softClipSortWindowSize, readPairSortWindowSize);
		// add safety margin to window size due the number of edge cases encountered which
		// result in a larger window size required
		return (int)(1.1 * windowSize);
	}
	public ReadPairConcordanceCalculator getReadPairConcordanceCalculator() {
		if (rpcc == null) {
			rpcc = ReadPairConcordanceCalculator.create(
					rpcMinFragmentSize == null ? 0 : rpcMinFragmentSize,
					rpcMaxFragmentSize == null ? 0 : rpcMaxFragmentSize,
					rpcConcordantPercentage,
					getMetrics().getInsertSizeDistribution(), getMetrics().getIdsvMetrics());
		}
		return rpcc;
	}
	@Override
	public int getMaxConcordantFragmentSize() {
		int fs = getMaxReadLength();
		fs = Math.max(getMaxReadMappedLength(), fs);
		if (getReadPairConcordanceCalculator() != null) {
			fs = Math.max(getReadPairConcordanceCalculator().maxConcordantFragmentSize(), fs);
		}
		return fs;
	}
	@Override
	public int getMinConcordantFragmentSize() {
		int fs = 1; // TODO: MIN_READ_LENGTH
		if (getReadPairConcordanceCalculator() != null) {
			fs = Math.max(getReadPairConcordanceCalculator().minConcordantFragmentSize(), fs);
		}
		return fs;
	}
	@Override
	public int getMaxReadLength() {
		return getMetrics().getIdsvMetrics().MAX_READ_LENGTH;
	}
	@Override
	public int getMaxReadMappedLength() {
		return getMetrics().getIdsvMetrics().MAX_READ_MAPPED_LENGTH;
	}
	public int getExpectedFragmentSize() {
		if (getMetrics().getInsertSizeMetrics() != null) {
			return Math.min((int)getMetrics().getInsertSizeMetrics().MEDIAN_INSERT_SIZE, getMaxConcordantFragmentSize());
		}
		return getMaxConcordantFragmentSize();
	}
	public GenomicProcessingContext getProcessContext() {
		return getContext();
	}
	public static CloseableIterator<DirectedEvidence> mergedIterator(List<SAMEvidenceSource> source, boolean parallel, EvidenceSortOrder eso) {
		List<CloseableIterator<DirectedEvidence>> toMerge = Lists.newArrayList();
		for (SAMEvidenceSource bam : source) {
			CloseableIterator<DirectedEvidence> it = bam.iterator(eso);
			if (parallel) {
				it = new AsyncBufferedIterator<>(it, bam.getFile() == null ? "" : bam.getFile().getName());
			}
			toMerge.add(it);
		}
		CloseableIterator<DirectedEvidence> merged = new AutoClosingMergedIterator<DirectedEvidence>(toMerge, eso == EvidenceSortOrder.EvidenceStartPosition ? DirectedEvidenceOrder.ByNatural : DirectedEvidenceOrder.BySAMStart);
		if (Defaults.SANITY_CHECK_DUMP_ITERATORS) {
			merged = new AutoClosingIterator<>(new DebugSpammingIterator<>(merged, "mergedIterator"));
		}
		return merged;
	}
	public static CloseableIterator<DirectedEvidence> mergedIterator(final List<SAMEvidenceSource> source, final QueryInterval[] intervals, EvidenceSortOrder eso) {
		List<CloseableIterator<DirectedEvidence>> toMerge = Lists.newArrayList();
		for (SAMEvidenceSource bam : source) {
			CloseableIterator<DirectedEvidence> it = bam.iterator(intervals, eso);
			toMerge.add(it);
		}
		CloseableIterator<DirectedEvidence> merged = new AutoClosingMergedIterator<DirectedEvidence>(toMerge,  eso == EvidenceSortOrder.EvidenceStartPosition ? DirectedEvidenceOrder.ByNatural : DirectedEvidenceOrder.BySAMStart);
		if (Defaults.SANITY_CHECK_DUMP_ITERATORS) {
			merged = new AutoClosingIterator<>(new DebugSpammingIterator<>(merged, "mergedIterator"));
		}
		return merged;
	}
	/**
	 * Maximum distance between the SAM alignment location of evidence, and the extrema of the
	 * breakend position supported by that evidence. 
	 * @return maximum out-of-order distance between evidence ordered by SAM alignment position and the breakend start position 
	 */
	public static int maximumWindowSize(ProcessingContext context, List<SAMEvidenceSource> sources, List<AssemblyEvidenceSource> assembly) {
		int maxSize = 0;
		for (EvidenceSource source : sources) {
			SAMEvidenceSource samSource = (SAMEvidenceSource)source;
			maxSize = Math.max(samSource.getMaxConcordantFragmentSize(), Math.max(samSource.getMaxReadLength(), samSource.getMaxReadMappedLength()));
		}
		if (assembly != null) {
			int maxAssemblyLength = assembly.stream().mapToInt(a -> a.getMaxAssemblyLength()).max().orElse(0);
			maxSize = Math.max(maxSize, maxAssemblyLength) + context.getVariantCallingParameters().maxBreakendHomologyLength;
		}
		return maxSize + 2 * (context.getVariantCallingParameters().breakendMargin + 1);
	}
	/**
	 * Provides hint as to whether the input file is known to contain only single-ended reads.
	 * @return
	 */
	public boolean knownSingleEnded() {
		return false;
	}
	public boolean isLongReadLibrary() {
		return getMetrics().getIdsvMetrics().MAX_READ_LENGTH >= getContext().getAssemblyParameters().longReadReadLengthThreshold;
	}
}