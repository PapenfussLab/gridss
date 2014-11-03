package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.concurrent.ForkJoinTask;

import au.edu.wehi.idsv.debruijn.anchored.DeBruijnAnchoredAssembler;
import au.edu.wehi.idsv.debruijn.subgraph.DeBruijnSubgraphAssembler;
import au.edu.wehi.idsv.pipeline.SortRealignedAssemblies;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.util.ParallelRunnableAction;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;


public class AssemblyEvidenceSource extends EvidenceSource {
	private static final Log log = Log.getInstance(AssemblyEvidenceSource.class);
	private final List<SAMEvidenceSource> source;
	private final int maxSourceFragSize;
	/**
	 * Generates assembly evidence based on the given evidence
	 * @param evidence evidence for creating assembly
	 * @param intermediateFileLocation location to store intermediate files
	 */
	public AssemblyEvidenceSource(ProcessingContext processContext, List<SAMEvidenceSource> evidence, File intermediateFileLocation) {
		super(processContext, intermediateFileLocation);
		this.source = evidence;
		int max = 0;
		for (SAMEvidenceSource s : evidence) {
			max = Math.max(max, s.getMaxConcordantFragmentSize());
		}
		maxSourceFragSize = max;
	}
	public void ensureAssembled() {
		if (!isProcessingComplete()) {
			process();
		}
		if (isRealignmentComplete()) {
			SortRealignedAssemblies step = new SortRealignedAssemblies(processContext, this);
			step.process(EnumSet.of(ProcessStep.SORT_REALIGNED_ASSEMBLIES));
			step.close();
		}
	}
	@Override
	protected Iterator<DirectedEvidence> perChrIterator(String chr) {
		FileSystemContext fsc = processContext.getFileSystemContext();
		@SuppressWarnings({ "unchecked", "rawtypes" }) // is there a less dirty way to do this downcast?
		Iterator<DirectedEvidence> downcast = (Iterator<DirectedEvidence>)(Iterator)iterator( 
				fsc.getAssemblyVcfForChr(input, chr),
				fsc.getRealignmentBamForChr(input, chr));
		return downcast;
	}
	@Override
	protected Iterator<DirectedEvidence> singleFileIterator() {
		FileSystemContext fsc = processContext.getFileSystemContext();
		@SuppressWarnings({ "unchecked", "rawtypes" })
		Iterator<DirectedEvidence> downcast = (Iterator<DirectedEvidence>)(Iterator)iterator(
			fsc.getAssemblyVcf(input),
			fsc.getRealignmentBam(input));
		return downcast;
	}
	private Iterator<VariantContextDirectedEvidence> iterator(File vcf, File realignment) {
		if (!isProcessingComplete()) {
			log.error("Assemblies not yet generated.");
			throw new RuntimeException("Assemblies not yet generated");
		}
		Iterator<SAMRecord> realignedIt; 
		if (isRealignmentComplete()) {
			SamReader realignedReader = processContext.getSamReader(realignment);
			realignedIt = new AutoClosingIterator<SAMRecord>(processContext.getSamReaderIterator(realignedReader), ImmutableList.<Closeable>of(realignedReader));
		} else {
			log.debug(String.format("Assembly realignment for %s not completed", vcf));
			realignedIt = ImmutableList.<SAMRecord>of().iterator();
		}
		VCFFileReader reader = new VCFFileReader(vcf, false);
		Iterator<VariantContextDirectedEvidence> it = new VariantContextDirectedEvidenceIterator(
				processContext,
				this,
				new AutoClosingIterator<VariantContext>(reader.iterator(), ImmutableList.<Closeable>of(reader)),
				realignedIt);
		// Change sort order from VCF sorted order to evidence position order
		it = new DirectEvidenceWindowedSortingIterator<VariantContextDirectedEvidence>(
				processContext,
				(int)((2 + processContext.getAssemblyParameters().maxSubgraphFragmentWidth + processContext.getAssemblyParameters().subgraphAssemblyMargin) * maxSourceFragSize),
				it);
		// FIXME: TODO: add remote assemblies to iterator
		return it;
	}
	private boolean isProcessingComplete() {
		boolean done = true;
		FileSystemContext fsc = processContext.getFileSystemContext();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				done &= IntermediateFileUtil.checkIntermediate(fsc.getAssemblyVcfForChr(input, seq.getSequenceName()));
				if (!done) return false;
				done &= IntermediateFileUtil.checkIntermediate(fsc.getRealignmentFastqForChr(input, seq.getSequenceName()));
				if (!done) return false;
			}
		} else {
			done &= IntermediateFileUtil.checkIntermediate(fsc.getAssemblyVcf(input));
			if (!done) return false;
			done &= IntermediateFileUtil.checkIntermediate(fsc.getRealignmentFastq(input));
			if (!done) return false;
		}
		return done;
	}
	protected void process() {
		log.info("START evidence assembly ", input);
		final Collection<Closeable> toClose = Collections.synchronizedCollection(new ArrayList<Closeable>());
		final SAMSequenceDictionary dict = processContext.getReference().getSequenceDictionary();
		try {
			if (processContext.shouldProcessPerChromosome()) {
				final List<Runnable> workers = Lists.newArrayList();
				for (int i = 0; i < dict.size(); i++) {
					final String seq = dict.getSequence(i).getSequenceName();
					workers.add(new Runnable() {
						@Override
						public void run() {
							List<Iterator<DirectedEvidence>> toMerge = Lists.newArrayList();
							for (SAMEvidenceSource bam : source) {
								Iterator<DirectedEvidence> it = bam.iterator(seq);
								if (it instanceof Closeable) toClose.add((Closeable)it);
								toMerge.add(it);
							}
							Iterator<DirectedEvidence> merged = Iterators.mergeSorted(toMerge, DirectedEvidenceOrder.ByNatural);
							new ContigAssembler(merged, processContext.getFileSystemContext().getAssemblyVcfForChr(input, seq), processContext.getFileSystemContext().getRealignmentFastqForChr(input, seq)).run();
							for (Iterator<DirectedEvidence> x : toMerge) {
								CloserUtil.close(x);
							}
						}
					});
				}
				if (ForkJoinTask.inForkJoinPool()) {
					log.info("Performing multithreaded assembly");
					ForkJoinTask.invokeAll(new ParallelRunnableAction(workers));
				} else {
					log.info("Performing singlethreaded assembly");
					for (Runnable r : workers) {
						r.run();
					}
				}
			} else {
				List<Iterator<DirectedEvidence>> toMerge = Lists.newArrayList();
				for (SAMEvidenceSource bam : source) {
					Iterator<DirectedEvidence> it = bam.iterator();
					if (it instanceof Closeable) toClose.add((Closeable)it);
					toMerge.add(it);
				}
				Iterator<DirectedEvidence> merged = Iterators.mergeSorted(toMerge, DirectedEvidenceOrder.ByNatural);
				new ContigAssembler(merged, processContext.getFileSystemContext().getAssemblyVcf(input), processContext.getFileSystemContext().getRealignmentFastq(input)).run();
				for (Iterator<DirectedEvidence> x : toMerge) {
					CloserUtil.close(x);
				}
			}
			log.info("SUCCESS evidence assembly ", input);
		} finally {
			for (Closeable x : toClose) {
				try {
					x.close();
				} catch (IOException e) {
					log.info(e, "Error closing stream");
				}
			}
		}
	}
	private class ContigAssembler implements Runnable {
		private Iterator<DirectedEvidence> it;
		private File breakendVcf;
		private File realignmentFastq;
		private FastqBreakpointWriter fastqWriter = null;
		private VariantContextWriter vcfWriter = null;
		private Queue<VariantContextDirectedEvidence> resortBuffer = new PriorityQueue<VariantContextDirectedEvidence>(32, IdsvVariantContext.ByLocationStart);
		private long maxAssembledPosition = Long.MIN_VALUE;
		private long lastFlushedPosition = Long.MIN_VALUE;
		private long lastProgress = 0;
		public ContigAssembler(Iterator<DirectedEvidence> it, File breakendVcf, File realignmentFastq) {
			this.it = it;
			this.breakendVcf = breakendVcf;
			this.realignmentFastq = realignmentFastq;
		}
		@Override
		public void run() {
			try {
				vcfWriter = processContext.getVariantContextWriter(breakendVcf);
				fastqWriter = new FastqBreakpointWriter(processContext.getFastqWriterFactory().newWriter(realignmentFastq));
				ReadEvidenceAssembler assembler = getAssembler();
				final ProgressLogger progress = new ProgressLogger(log);
				while (it.hasNext()) {
					DirectedEvidence readEvidence = it.next();
					if (readEvidence instanceof NonReferenceReadPair) {
						progress.record(((NonReferenceReadPair)readEvidence).getLocalledMappedRead());
					} else if (readEvidence instanceof SoftClipEvidence) {
						progress.record(((SoftClipEvidence)readEvidence).getSAMRecord());
					}
					// Need to process assembly evidence first since assembly calls are made when the
					// evidence falls out of scope so processing a given position will emit evidence
					// for a previous position (for it is only at this point we know there is no more
					// evidence for the previous position).
					processAssembly(assembler.addEvidence(readEvidence), resortBuffer, fastqWriter, vcfWriter);
					
					if (maxAssembledPosition / 1000000 > lastProgress / 1000000) {
						lastProgress = maxAssembledPosition;
						log.info(String.format("Assembly at %s:%d %s",
								processContext.getDictionary().getSequence(processContext.getLinear().getReferenceIndex(lastProgress)).getSequenceName(),
								processContext.getLinear().getReferencePosition(lastProgress),
								assembler.getStateSummaryMetrics()));
					}
				}
				processAssembly(assembler.endOfEvidence(), resortBuffer, fastqWriter, vcfWriter);
				flushWriterQueueBefore(Long.MAX_VALUE, resortBuffer, fastqWriter, vcfWriter);
			} finally {
				if (fastqWriter != null) fastqWriter.close();
				if (vcfWriter != null) vcfWriter.close();
			}
		}
		private void processAssembly(
				Iterable<VariantContextDirectedEvidence> evidenceList,
				Queue<VariantContextDirectedEvidence> buffer,
				FastqBreakpointWriter fastqWriter,
				VariantContextWriter vcfWriter) {
	    	if (evidenceList != null) {
		    	for (VariantContextDirectedEvidence a : evidenceList) {
		    		if (a != null) {
			    		maxAssembledPosition = Math.max(maxAssembledPosition, processContext.getLinear().getLinearCoordinate(a.getReferenceIndex(), a.getStart()));
			    		if (Defaults.WRITE_FILTERED_CALLS || !a.isFiltered()) {
			    			buffer.add(a);
			    		}
		    		}
		    	}
	    	}
	    	flushWriterQueueBefore(
	    			maxAssembledPosition - (long)((processContext.getAssemblyParameters().maxSubgraphFragmentWidth * 2) * getMaxConcordantFragmentSize()),
	    			buffer,
	    			fastqWriter,
	    			vcfWriter);
	    }
		private void flushWriterQueueBefore(
				long flushBefore,
				Queue<VariantContextDirectedEvidence> buffer,
				FastqBreakpointWriter fastqWriter,
				VariantContextWriter vcfWriter) {
			while (!buffer.isEmpty() && processContext.getLinear().getLinearCoordinate(buffer.peek().getReferenceIndex(), buffer.peek().getStart()) < flushBefore) {
				long pos = processContext.getLinear().getLinearCoordinate(buffer.peek().getReferenceIndex(), buffer.peek().getStart());
				VariantContextDirectedEvidence evidence = buffer.poll();
				if (pos < lastFlushedPosition) {
					log.error(String.format("Sanity check failure: assembly breakend %s written out of order.", evidence.getID()));
				}
				lastFlushedPosition = pos;
				if (Defaults.WRITE_FILTERED_CALLS || !evidence.isFiltered()) {
					vcfWriter.add(evidence);
				}
				if (!evidence.isFiltered() && processContext.getRealignmentParameters().shouldRealignBreakend(evidence)) {
					fastqWriter.write(evidence);
				}
			}
		}
	}
	private ReadEvidenceAssembler getAssembler() {
    	switch (processContext.getAssemblyParameters().method) {
	    	case DEBRUIJN_PER_POSITION:
	    		return new DeBruijnAnchoredAssembler(processContext, this);
	    	case DEBRUIJN_SUBGRAPH:
	    		return new DeBruijnSubgraphAssembler(processContext, this);
	    	default:
	    		throw new IllegalArgumentException("Unknown assembly method.");
    	}
    }
	@Override
	public int getMaxConcordantFragmentSize() {
		return maxSourceFragSize;
	}
}
