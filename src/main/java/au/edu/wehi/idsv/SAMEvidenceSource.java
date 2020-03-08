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
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import gridss.ComputeSamTags;
import gridss.ExtractSVReads;
import gridss.SoftClipsToSplitReads;
import gridss.analysis.CollectGridssMetrics;
import gridss.analysis.StructuralVariantReadMetrics;
import gridss.cmdline.CommandLineProgramHelper;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;

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
	private StructuralVariantReadMetrics svMetrics;
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
	public StructuralVariantReadMetrics getSVMetrics() {
		File svMetricsFile = getContext().getFileSystemContext().getSVMetrics(getFile());
		if (svMetrics == null && svMetricsFile.exists()) {
			for (StructuralVariantReadMetrics metric : Iterables.filter(MetricsFile.readBeans(svMetricsFile), StructuralVariantReadMetrics.class)) {
				if (metric.SAMPLE == null && metric.LIBRARY == null && metric.READ_GROUP == null) {
					svMetrics = metric;
					break;
				}
			}
		}
		return svMetrics;
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
				List<String> args = Lists.newArrayList(
						"INPUT=" + getFile().getPath(),
						"OUTPUT=" + getContext().getFileSystemContext().getMetricsPrefix(getFile()).getPath(),
						"THRESHOLD_COVERAGE=" + getContext().getConfig().maxCoverage,
						"FILE_EXTENSION=null",
						"GRIDSS_PROGRAM=null",
						"GRIDSS_PROGRAM=CollectCigarMetrics",
						"GRIDSS_PROGRAM=CollectMapqMetrics",
						"GRIDSS_PROGRAM=CollectTagMetrics",
						"GRIDSS_PROGRAM=CollectIdsvMetrics",
						"GRIDSS_PROGRAM=ReportThresholdCoverage",
						"PROGRAM=null");
				if (!knownSingleEnded()) {
					// Don't run CollectInsertSizeMetrics
					args.add("PROGRAM=CollectInsertSizeMetrics");
				} else {
					// The CollectMultipleMetrics super class complains if no PROGRAM set so
					// we'll just collect some stuff that is useful, but we don't actually
					// use yet
					args.add("PROGRAM=CollectAlignmentSummaryMetrics");
				}
				if (getContext().getCalculateMetricsRecordCount() < Integer.MAX_VALUE) {
					args.add("STOP_AFTER=" + getContext().getCalculateMetricsRecordCount());
				}
				execute(new CollectGridssMetrics(), args);
			}
			metrics = new IdsvSamFileMetrics(getContext(), getFile(), knownSingleEnded());
		}
	}
	protected void execute(CommandLineProgram cmd, List<String> args) {
		if (cmd instanceof ReferenceCommandLineProgram) {
			((ReferenceCommandLineProgram) cmd).setReference(getContext().getReference());
		}
		if (getContext().getCommandLineProgram() == null) {
			args.add("REFERENCE_SEQUENCE=" + getContext().getReferenceFile());
			args.add("TMP_DIR=" + getContext().getFileSystemContext().getTemporaryDirectory());
			args.add("MAX_RECORDS_IN_RAM=" + getContext().getFileSystemContext().getMaxBufferedRecordsPerFile());
		} else {
			args.addAll(CommandLineProgramHelper.getCommonArgs(getContext().getCommandLineProgram()));
		}
		int result = cmd.instanceMain(args.toArray(new String[] {}));
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
							List<String> args = Lists.newArrayList(
									"INPUT=" + in.getPath(),
									"OUTPUT=" + extractedFile.getPath(),
									"UNMAPPED_READS=false", // saves intermediate file space
									"METRICS_OUTPUT=" + getContext().getFileSystemContext().getSVMetrics(getFile()),
									"MIN_CLIP_LENGTH=" + getContext().getConfig().getSoftClip().minLength,
									"INSERT_SIZE_METRICS=" + getContext().getFileSystemContext().getInsertSizeMetrics(getFile()),
									// Picard tools does not mark duplicates correctly. We need to keep them so we can
									// fix the duplicate marking in ComputeSamTags
									"INCLUDE_DUPLICATES=true");
							if (rpcMinFragmentSize != null) args.add("READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE=" + rpcMinFragmentSize);
							if (rpcMaxFragmentSize != null) args.add("READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE=" + rpcMaxFragmentSize);
							if (rpcConcordantPercentage != null) args.add("READ_PAIR_CONCORDANT_PERCENT=" + rpcConcordantPercentage);
							execute(new ExtractSVReads(), args);
						}
						SAMFileUtil.sort(getContext().getFileSystemContext(), extractedFile, querysortedFile, SortOrder.queryname);
						if (gridss.Defaults.DELETE_TEMPORARY_FILES) {
							FileHelper.delete(extractedFile, true);
						}
					}
					log.info("Computing SAM tags for " + svFile);
					List<String> args = Lists.newArrayList(
							"INPUT=" + querysortedFile.getPath(),
							"OUTPUT=" + taggedFile.getPath());
					execute(new ComputeSamTags(), args);
					if (gridss.Defaults.DELETE_TEMPORARY_FILES) {
						FileHelper.delete(querysortedFile, true);
					}
				}
				log.info("Identifying split reads for " + getFile().getAbsolutePath());
				List<String> args = Lists.newArrayList(
						"WORKER_THREADS=" + getProcessContext().getWorkerThreadCount(),
						"INPUT=" + taggedFile.getPath(),
						"OUTPUT=" + withsplitreadsFile.getPath(),
						"REALIGN_EXISTING_SPLIT_READS=" + Boolean.toString(getContext().getConfig().getSoftClip().realignSplitReads));
						// realignment.* not soft-clip
						//"MIN_CLIP_LENGTH=" + getContext().getConfig().
						//"MIN_CLIP_QUAL=" + getContext().getConfig().getSoftClip().minAverageQual);
				SoftClipsToSplitReads program =  new SoftClipsToSplitReads();
				program.setReference(getProcessContext().getReference());
				program.setFileSystemContext(getProcessContext().getFileSystemContext());
				execute(program, args);
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
		SAMRecordIterator it = tryOpenReader(reader, queryInterval.asQueryInterval());
		Iterator<DirectedEvidence> eit = asEvidence(it, eso);
		eit = Iterators.filter(eit, e -> QueryIntervalUtil.overlaps(intervals, e.getBreakendSummary()));
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
		it = Iterators.transform(it, r -> transform(r));
		it = Iterators.filter(it, r -> !shouldFilter(r));		
		Iterator<DirectedEvidence> eit = new DirectedEvidenceIterator(it, this, minIndelSize());
		eit = Iterators.filter(eit, e -> !shouldFilter(e));
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
		if (SAMRecordUtil.isDovetailing(r, PairOrientation.FR, getContext().getConfig().dovetailMargin)) {
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
				return getContext().getBlacklistedRegions();
			}
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
			if (average(sce.getBreakendQuality()) < scc.minAverageQual) return true;
			if (config.adapters.isAdapterSoftClip(sce)) return true;
			if (!AssemblyAttributes.isAssembly(sce.getSAMRecord())) {
				// TODO: symmetrical identity and entropy filters on both sides
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
			rpcc = ReadPairConcordanceCalculator.create(rpcMinFragmentSize, rpcMaxFragmentSize, rpcConcordantPercentage, getMetrics().getInsertSizeDistribution(), getMetrics().getIdsvMetrics());
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
		return merged;
	}
	public static CloseableIterator<DirectedEvidence> mergedIterator(final List<SAMEvidenceSource> source, final QueryInterval[] intervals, EvidenceSortOrder eso) {
		List<CloseableIterator<DirectedEvidence>> toMerge = Lists.newArrayList();
		for (SAMEvidenceSource bam : source) {
			CloseableIterator<DirectedEvidence> it = bam.iterator(intervals, eso);
			toMerge.add(it);
		}
		CloseableIterator<DirectedEvidence> merged = new AutoClosingMergedIterator<DirectedEvidence>(toMerge,  eso == EvidenceSortOrder.EvidenceStartPosition ? DirectedEvidenceOrder.ByNatural : DirectedEvidenceOrder.BySAMStart);
		return merged;
	}
	/**
	 * Maximum distance between the SAM alignment location of evidence, and the extrema of the
	 * breakend position supported by that evidence. 
	 * @return maximum out-of-order distance between evidence ordered by SAM alignment position and the breakend start position 
	 */
	public static int maximumWindowSize(ProcessingContext context, List<SAMEvidenceSource> sources, AssemblyEvidenceSource assembly) {
		int maxSize = 0;
		for (EvidenceSource source : sources) {
			SAMEvidenceSource samSource = (SAMEvidenceSource)source;
			maxSize = Math.max(samSource.getMaxConcordantFragmentSize(), Math.max(samSource.getMaxReadLength(), samSource.getMaxReadMappedLength()));
		}
		if (assembly != null) {
			maxSize = Math.max(maxSize, assembly.getMaxAssemblyLength()) + context.getVariantCallingParameters().maxBreakendHomologyLength;
		}
		return maxSize + 2 * (context.getVariantCallingParameters().breakendMargin + 1);
	}
	/**
	 * Provides hint as to whether the input file is known to contain only single-ended reads.
	 * @return
	 */
	public boolean knownSingleEnded() { return false; }
}