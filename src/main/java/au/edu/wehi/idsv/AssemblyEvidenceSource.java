package au.edu.wehi.idsv;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ExecutorService;

import com.google.common.collect.Iterators;

import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.configuration.AssemblyConfiguration;
import au.edu.wehi.idsv.debruijn.positional.PositionalAssembler;
import au.edu.wehi.idsv.sam.SAMFileUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.FileHelper;
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
		super(processContext, assemblyFile, -1);
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
		SAMFileHeader header = getContext().getBasicSamHeader();
		header.setSortOrder(SortOrder.unsorted);
		// TODO: add assembly @PG header
		File tmpout = FileSystemContext.getWorkingFileFor(getFile());
		File filteredout = FileSystemContext.getWorkingFileFor(getFile(), "filtered.");
		try (SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, tmpout)) {
			try (SAMFileWriter filteredWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, filteredout)) {
				Iterator<SAMRecord> it = getAllAssemblies(threadpool);
				while (it.hasNext()) {
					SAMRecord asm = it.next();
					// some assemblies actually match the reference and we can ignore these
					// such assemblies are caused by sequencing errors or SNVs causing
					// spurious soft clips
					SAMRecordUtil.unclipExactReferenceMatches(getContext().getReference(), asm);
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
			if (!getContext().getAssemblyParameters().writeFiltered) {
				FileHelper.delete(filteredout, true);
			}
		}
		if (filteredout.exists()) {
			SAMFileUtil.sort(getContext().getFileSystemContext(), filteredout, filteredout, SortOrder.coordinate);
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
	public boolean shouldFilter(SAMRecord r) {
		return shouldFilterAssembly(r) || super.shouldFilter(r);
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
	private Iterator<SAMRecord> getAllAssemblies(ExecutorService threadpool) {
		// TODO multi-threaded version
		// Chunk into tasks
		// Schedule tasks on thread pool
		// Get records from blocking output queue
			// stop when no outstanding assembly tasks
		if (threadpool == null) {
			return getAllAssemblies_single_threaded();
		}
		return getAllAssemblies_single_threaded();
	}
	private Iterator<SAMRecord> getAllAssemblies_single_threaded() {
		ProgressLogger progressLog = new ProgressLogger(log);
		List<Iterator<SAMRecord>> list = new ArrayList<>();
		for (BreakendDirection direction : BreakendDirection.values()) {
			CloseableIterator<DirectedEvidence> it = mergedIterator(source);
			Iterator<DirectedEvidence> throttledIt = throttled(it);
			ProgressLoggingDirectedEvidenceIterator<DirectedEvidence> loggedIt = new ProgressLoggingDirectedEvidenceIterator<>(getContext(), throttledIt, progressLog);
			Iterator<SAMRecord> evidenceIt = new PositionalAssembler(getContext(), this, loggedIt, direction);
	    	list.add(evidenceIt);
		}
		return Iterators.concat(list.iterator());
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
}
