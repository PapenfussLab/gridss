package au.edu.wehi.idsv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import java.io.File;
import java.util.Iterator;
import java.util.List;

import au.edu.wehi.idsv.debruijn.anchored.DeBruijnAnchoredAssembler;
import au.edu.wehi.idsv.debruijn.subgraph.DeBruijnSubgraphAssembler;
import au.edu.wehi.idsv.metrics.CompositeMetrics;
import au.edu.wehi.idsv.metrics.RelevantMetrics;

import com.google.common.base.Function;
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;


public class AssemblyEvidenceSource extends EvidenceSource {
	private static final Log log = Log.getInstance(AssemblyEvidenceSource.class);
	private final List<SAMEvidenceSource> source;
	private final RelevantMetrics metrics;
	/**
	 * Generates assembly evidence based on the given evidence
	 * @param evidence evidence for creating assembly
	 * @param intermediateFileLocation location to store intermediate files
	 */
	public AssemblyEvidenceSource(ProcessingContext processContext, List<SAMEvidenceSource> evidence, File intermediateFileLocation) {
		super(processContext, intermediateFileLocation);
		this.source = evidence;
		this.metrics = new CompositeMetrics(Iterables.transform(evidence, new Function<SAMEvidenceSource, RelevantMetrics>() {
			@Override
			public RelevantMetrics apply(SAMEvidenceSource arg) {
				return arg.getMetrics();
			}
		}));
	}
	public void ensureAssembled() {
		if (!isProcessingComplete()) {
			log.info("START evidence assembly ", input);
			process();
			log.info("SUCCESS evidence assembly ", input);
		}
	}
	@Override
	public RelevantMetrics getMetrics() {
		return metrics;
	}
	@Override
	protected DirectedEvidenceFileIterator perChrIterator(String chr) {
		FileSystemContext fsc = processContext.getFileSystemContext();
		return new DirectedEvidenceFileIterator(processContext, this,
				null,
				null,
				isRealignmentComplete() ? fsc.getRealignmentBamForChr(input, chr) : null,
				fsc.getBreakendVcfForChr(input, chr));
	}
	@Override
	protected DirectedEvidenceFileIterator singleFileIterator() {
		FileSystemContext fsc = processContext.getFileSystemContext();
		return new DirectedEvidenceFileIterator(processContext, this,
				null,
				null,
				isRealignmentComplete() ? fsc.getRealignmentBam(input) : null,
				fsc.getBreakendVcf(input));
	}
	@Override
	protected boolean isProcessingComplete() {
		boolean done = true;
		FileSystemContext fsc = processContext.getFileSystemContext();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				done &= IntermediateFileUtil.checkIntermediate(fsc.getBreakendVcfForChr(input, seq.getSequenceName()));
				done &= IntermediateFileUtil.checkIntermediate(fsc.getRealignmentFastqForChr(input, seq.getSequenceName()));
			}
		} else {
			done &= IntermediateFileUtil.checkIntermediate(fsc.getBreakendVcf(input));
			done &= IntermediateFileUtil.checkIntermediate(fsc.getRealignmentFastq(input));
		}
		return done;
	}
	@Override
	protected void process() {
		List<CloseableIterator<DirectedEvidence>> toClose = Lists.newArrayList();
		SAMSequenceDictionary dict = processContext.getReference().getSequenceDictionary();
		try {
			if (processContext.shouldProcessPerChromosome()) {
				for (int i = 0; i < dict.size(); i++) {
					String seq = dict.getSequence(i).getSequenceName();
					List<CloseableIterator<DirectedEvidence>> toMerge = Lists.newArrayList();
					for (SAMEvidenceSource bam : source) {
						CloseableIterator<DirectedEvidence> it = bam.iterator(seq);
						toClose.add(it);
						toMerge.add(it);
					}
					Iterator<DirectedEvidence> merged = Iterators.mergeSorted(toMerge, DirectedEvidenceOrder.ByNatural);
					assemble(merged, processContext.getFileSystemContext().getBreakendVcfForChr(input, seq), processContext.getFileSystemContext().getRealignmentFastqForChr(input, seq));
					for (CloseableIterator<DirectedEvidence> x : toMerge) {
						x.close();
					}
				}
			} else {
				List<CloseableIterator<DirectedEvidence>> toMerge = Lists.newArrayList();
				for (SAMEvidenceSource bam : source) {
					CloseableIterator<DirectedEvidence> it = bam.iterator();
					toClose.add(it);
					toMerge.add(it);
				}
				Iterator<DirectedEvidence> merged = Iterators.mergeSorted(toMerge, DirectedEvidenceOrder.ByNatural);
				assemble(merged, processContext.getFileSystemContext().getBreakendVcf(input), processContext.getFileSystemContext().getRealignmentFastq(input));
				for (CloseableIterator<DirectedEvidence> x : toMerge) {
					x.close();
				}
			}
		} finally {
			for (CloseableIterator<DirectedEvidence> x : toClose) {
				x.close();
			}
		}
	}
	private void assemble(Iterator<DirectedEvidence> it, File breakendVcf, File realignmentFastq) {
		FastqBreakpointWriter fastqWriter = null;
		VariantContextWriter vcfWriter = null;
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
				processAssembly(assembler.addEvidence(readEvidence), fastqWriter, vcfWriter);
			}
			processAssembly(assembler.endOfEvidence(), fastqWriter, vcfWriter);
		} finally {
			if (fastqWriter != null) fastqWriter.close();
			if (vcfWriter != null) vcfWriter.close();
		}
	}
	private void processAssembly(Iterable<VariantContextDirectedBreakpoint> evidenceList, FastqBreakpointWriter fastqWriter, VariantContextWriter vcfWriter) {
    	if (evidenceList != null) {
	    	for (VariantContextDirectedBreakpoint a : evidenceList) {
	    		processAssembly(a, fastqWriter, vcfWriter);
	    	}
    	}
    }
	private void processAssembly(VariantContextDirectedBreakpoint evidence, FastqBreakpointWriter fastqWriter, VariantContextWriter vcfWriter) {
		if (evidence != null) {
			vcfWriter.add(evidence); // write out failed assemblies as well for debugging purposes
			if (!evidence.isFiltered() && processContext.getRealignmentParameters().shouldRealignBreakend(evidence)) {
				fastqWriter.write(evidence);
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
}
