package au.edu.wehi.idsv;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.configuration.GridssConfiguration;
import au.edu.wehi.idsv.configuration.SoftClipConfiguration;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.sam.SAMFileUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.util.AutoClosingMergedIterator;
import au.edu.wehi.idsv.util.BufferedIterator;
import au.edu.wehi.idsv.util.FileHelper;
import au.edu.wehi.idsv.validation.OrderAssertingIterator;
import gridss.ComputeSamTags;
import gridss.ExtractSVReads;
import gridss.SoftClipsToSplitReads;
import gridss.analysis.CollectGridssMetrics;
import gridss.analysis.StructuralVariantReadMetrics;
import gridss.cmdline.CommandLineProgramHelper;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;

/**
 * Structural variation evidence based on read pairs from a single SAM/BAM.  
 * @author Daniel Cameron
 *
 */
public class SAMEvidenceSource extends EvidenceSource {
	private static final Log log = Log.getInstance(SAMEvidenceSource.class);
	protected final SamReaderFactory factory = SamReaderFactory.makeDefault();
	private final int sourceCategory;
	private final ReadPairConcordanceMethod rpcMethod;
	private final int rpcMinFragmentSize;
	private final int rpcMaxFragmentSize;
	private final double rpcConcordantPercentage;
	private IdsvSamFileMetrics metrics;
	private StructuralVariantReadMetrics svMetrics;
	private ReadPairConcordanceCalculator rpcc;
	public SAMEvidenceSource(ProcessingContext processContext, File file, File nameSorted, int sourceCategory) {
		this(processContext, file, nameSorted, sourceCategory, ReadPairConcordanceMethod.SAM_FLAG, 0, 0, 0);
	}
	public SAMEvidenceSource(ProcessingContext processContext, File file, File nameSorted, int sourceCategory, int minFragmentSize, int maxFragmentSize) {
		this(processContext, file, nameSorted, sourceCategory, ReadPairConcordanceMethod.FIXED, minFragmentSize, maxFragmentSize, 0);
	}
	public SAMEvidenceSource(ProcessingContext processContext, File file, File nameSorted, int sourceCategory, double concordantPercentage) {
		this(processContext, file, nameSorted, sourceCategory, ReadPairConcordanceMethod.PERCENTAGE, 0, 0, concordantPercentage);
	}
	protected SAMEvidenceSource(ProcessingContext processContext, File file, File nameSorted,
			int sourceCategory, ReadPairConcordanceMethod rpcMethod, int rpcMinFragmentSize, int rpcMaxFragmentSize, double rpcConcordantPercentage) {
		super(processContext, file, nameSorted);
		this.sourceCategory = sourceCategory;
		this.rpcMethod = rpcMethod;
		this.rpcMinFragmentSize = rpcMinFragmentSize;
		this.rpcMaxFragmentSize = rpcMaxFragmentSize;
		this.rpcConcordantPercentage = rpcConcordantPercentage;
	}
	public IdsvSamFileMetrics getMetrics() {
		ensureMetrics();
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
	public void ensureMetrics() {
		if (metrics == null) {
			File idsvFile = getContext().getFileSystemContext().getIdsvMetrics(getFile());
			File cigarFile = getContext().getFileSystemContext().getCigarMetrics(getFile());
			File mapqFile = getContext().getFileSystemContext().getMapqMetrics(getFile());
			if (!idsvFile.exists() || !cigarFile.exists() || !mapqFile.exists()) {
				log.info("Calculating metrics for " + getFile().getAbsolutePath());
				List<String> args = Lists.newArrayList(
						"INPUT=" + getFile().getAbsolutePath(),
						"OUTPUT=" + getContext().getFileSystemContext().getMetricsPrefix(getFile()).getAbsolutePath(),
						"THRESHOLD_COVERAGE=" + getContext().getConfig().maxCoverage,
						"FILE_EXTENSION=null",
						"GRIDSS_PROGRAM=null",
						"GRIDSS_PROGRAM=CollectCigarMetrics",
						"GRIDSS_PROGRAM=CollectMapqMetrics",
						"GRIDSS_PROGRAM=CollectTagMetrics",
						"GRIDSS_PROGRAM=CollectIdsvMetrics",
						"GRIDSS_PROGRAM=ReportThresholdCoverage",
						// The CollectMultipleMetrics super class complains if no PROGRAM set so
						// we'll just collect some stuff that is useful, but we don't actually
						// use yet
						"PROGRAM=null",
						"PROGRAM=CollectAlignmentSummaryMetrics",
						"PROGRAM=QualityScoreDistribution");
				if (!knownSingleEnded()) {
					// Don't run CollectInsertSizeMetrics
					args.add("PROGRAM=CollectInsertSizeMetrics");
				}
				if (getContext().getCalculateMetricsRecordCount() < Integer.MAX_VALUE) {
					args.add("STOP_AFTER=" + getContext().getCalculateMetricsRecordCount());
				}
				execute(new CollectGridssMetrics(), args);
			}
			metrics = new IdsvSamFileMetrics(getContext(), getFile());
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
	public void ensureExtracted() throws IOException {
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
								if (getContext().getConfig().multimapping) {
									throw new IllegalArgumentException(String.format("Missing INPUT_NAME_SORTED for %s."
											+ " Both coordinate and name sorted input files must be supplied when multi-mapping mode is enabled.", getFile()));
								}
								in = getFile();
							}
							List<String> args = Lists.newArrayList(
									"INPUT=" + in.getAbsolutePath(),
									"OUTPUT=" + extractedFile.getAbsolutePath(),
									"UNMAPPED_READS=false", // saves intermediate file space
									"METRICS_OUTPUT=" + getContext().getFileSystemContext().getSVMetrics(getFile()),
									"MIN_CLIP_LENGTH=" + getContext().getConfig().getSoftClip().minLength,
									"READ_PAIR_CONCORDANCE_METHOD=" + rpcMethod.name(),
									"FIXED_READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE=" + rpcMinFragmentSize,
									"FIXED_READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE=" + rpcMaxFragmentSize,
									"READ_PAIR_CONCORDANT_PERCENT=" + rpcConcordantPercentage,
									"INSERT_SIZE_METRICS=" + getContext().getFileSystemContext().getInsertSizeMetrics(getFile()));
							execute(new ExtractSVReads(), args);
						}
						SAMFileUtil.sort(getContext().getFileSystemContext(), extractedFile, querysortedFile, SortOrder.queryname);
						if (gridss.Defaults.DELETE_TEMPORARY_FILES) {
							FileHelper.delete(extractedFile, true);
						}
					}
					log.info("Computing SAM tags for " + svFile);
					List<String> args = Lists.newArrayList(
							"INPUT=" + querysortedFile.getAbsolutePath(),
							"OUTPUT=" + taggedFile.getAbsolutePath());
					execute(new ComputeSamTags(), args);
					if (gridss.Defaults.DELETE_TEMPORARY_FILES) {
						FileHelper.delete(querysortedFile, true);
					}
				}
				log.info("Identifying split reads for " + getFile().getAbsolutePath());
				List<String> args = Lists.newArrayList(
						"WORKER_THREADS=" + getProcessContext().getWorkerThreadCount(),
						"INPUT=" + taggedFile.getAbsolutePath(),
						"OUTPUT=" + withsplitreadsFile.getAbsolutePath());
						// realignment.* not soft-clip
						//"MIN_CLIP_LENGTH=" + getContext().getConfig().
						//"MIN_CLIP_QUAL=" + getContext().getConfig().getSoftClip().minAverageQual);
				execute(new SoftClipsToSplitReads(), args);
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
	public CloseableIterator<DirectedEvidence> iterator(final QueryInterval interval) {
		SamReader reader = getReader();
		// expand query bounds as the alignment for a discordant read pair could fall before or after the breakend interval we are extracting
		SAMRecordIterator it = reader.queryOverlapping(new QueryInterval[] {
			new QueryInterval(interval.referenceIndex, interval.start - getMaxConcordantFragmentSize() - 1, interval.end + getMaxConcordantFragmentSize() + 1)});
		Iterator<DirectedEvidence> eit = asEvidence(it);
		final BreakendSummary bsf = new BreakendSummary(interval.referenceIndex, BreakendDirection.Forward, interval.start, interval.start, interval.end);
		final BreakendSummary bsb = new BreakendSummary(interval.referenceIndex, BreakendDirection.Backward, interval.start, interval.start, interval.end);
		eit = Iterators.filter(eit, e -> bsf.overlaps(e.getBreakendSummary()) || bsb.overlaps(e.getBreakendSummary()));
		return new AutoClosingIterator<>(eit, reader, it);
	}
	public CloseableIterator<DirectedEvidence> iterator() {
		SamReader reader = getReader();
		SAMRecordIterator it = reader.iterator();
		it.assertSorted(SortOrder.coordinate);
		Iterator<DirectedEvidence> eit = asEvidence(it);
		return new AutoClosingIterator<>(eit, reader, it);
	}
	private SamReader getReader() {
		File svFile = getContext().getFileSystemContext().getSVBam(getFile());
		SamReader reader = factory.open(svFile.exists() ? svFile : getFile());
		return reader;
	}
	private Iterator<DirectedEvidence> asEvidence(Iterator<SAMRecord> it) {
		it = new BufferedIterator<>(it, 2); // TODO: remove when https://github.com/samtools/htsjdk/issues/760 is resolved 
		it = Iterators.transform(it, r -> transform(r));
		it = Iterators.filter(it, r -> !shouldFilter(r));		
		Iterator<DirectedEvidence> eit = new DirectedEvidenceIterator(it, this, minIndelSize());
		eit = Iterators.filter(eit, e -> !shouldFilter(e));
		eit = new DirectEvidenceWindowedSortingIterator<DirectedEvidence>(getContext(), getSortWindowSize(), eit);
		if (Defaults.SANITY_CHECK_ITERATORS) {
			// Can't enforce pairing as there may actually be duplicates depending on how multi-mapping alignment was performed
			//eit = new PairedEvidenceTracker<DirectedEvidence>(getFile().getName(), eit);
			eit = new OrderAssertingIterator<DirectedEvidence>(eit, DirectedEvidenceOrder.ByNatural);
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
			String mc = r.getStringAttribute(SAMTag.MC.name());
			if (mc != null) {
				Cigar mateCigar = TextCigarCodec.decode(mc);
				mateEnd += mateCigar.getReferenceLength() - 1;
			}
			if (getBlacklistedRegions().overlaps(mateRef, mateStart, mateEnd)) {
				r.setMateUnmappedFlag(true);
			}
		}
		if (r.getAttribute(SAMTag.SA.name()) != null) {
			SAMSequenceDictionary dict = getContext().getDictionary();
			r.setAttribute(SAMTag.SA.name(), ChimericAlignment.getChimericAlignments(r).stream()
					.filter(ca -> !getBlacklistedRegions().overlaps(
							dict.getSequence(ca.rname).getSequenceIndex(),
							ca.pos,
							ca.pos + ca.cigar.getReferenceLength() - 1))
					.map(ca -> ca.toString())
					.collect(Collectors.joining(";")));
		}
		return r;
	}
	public boolean shouldFilter(SAMRecord r) {
		if (r.getReadUnmappedFlag()) {
			return true;
		}
		if (getContext().isFilterDuplicates() && r.getDuplicateReadFlag()) {
			return true;
		}
		if (SAMRecordUtil.isDovetailing(r, PairOrientation.FR, getContext().getConfig().dovetailMargin)) {
			return true;
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
				blacklist = IntervalBed.merge(getContext().getDictionary(), getContext().getLinear(), ImmutableList.of(
						getContext().getBlacklistedRegions(),
						new IntervalBed(getContext().getDictionary(), getContext().getLinear(), coverageBlacklist)
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
		// MMMMMMMMMMMMMMMMMM>
		// ^--read length --^
		// ^--------max concordant fragment size-----^
		//                   |-----breakend call-----|
		//                   |----------breakend call--------|
		//                                                   <SSSSSSSSSM
		//                   ^--------max concordant fragment size-----^
		// ^ alignment start                                           ^ alignment start
		int readPairSortWindowSize = getMaxConcordantFragmentSize() + getMaxReadMappedLength() + 1;
		int windowSize = Math.max(softClipSortWindowSize, readPairSortWindowSize);
		return windowSize;
	}
	public ReadPairConcordanceCalculator getReadPairConcordanceCalculator() {
		if (rpcc == null) {
			rpcc = ReadPairConcordanceCalculator.create(rpcMethod, rpcMinFragmentSize, rpcMaxFragmentSize, rpcConcordantPercentage, getMetrics().getInsertSizeDistribution(), getMetrics().getIdsvMetrics());
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
	public static CloseableIterator<DirectedEvidence> mergedIterator(List<SAMEvidenceSource> source, boolean parallel) {
		List<CloseableIterator<DirectedEvidence>> toMerge = Lists.newArrayList();
		for (SAMEvidenceSource bam : source) {
			CloseableIterator<DirectedEvidence> it = bam.iterator();
			if (parallel) {
				it = new AsyncBufferedIterator<>(it, bam.getFile() == null ? "" : bam.getFile().getName());
			}
			toMerge.add(it);
		}
		CloseableIterator<DirectedEvidence> merged = new AutoClosingMergedIterator<DirectedEvidence>(toMerge, DirectedEvidenceOrder.ByNatural);
		return merged;
	}
	public static CloseableIterator<DirectedEvidence> mergedIterator(final List<SAMEvidenceSource> source, final QueryInterval intervals) {
		List<CloseableIterator<DirectedEvidence>> toMerge = Lists.newArrayList();
		for (SAMEvidenceSource bam : source) {
			CloseableIterator<DirectedEvidence> it = bam.iterator(intervals);
			toMerge.add(it);
		}
		CloseableIterator<DirectedEvidence> merged = new AutoClosingMergedIterator<DirectedEvidence>(toMerge, DirectedEvidenceOrder.ByNatural);
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