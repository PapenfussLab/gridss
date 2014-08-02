package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import au.edu.wehi.idsv.debruijn.anchored.DeBruijnAnchoredAssembler;
import au.edu.wehi.idsv.debruijn.subgraph.DeBruijnSubgraphAssembler;
import au.edu.wehi.idsv.metrics.CompositeMetrics;
import au.edu.wehi.idsv.metrics.RelevantMetrics;
import au.edu.wehi.idsv.util.AutoClosingIterator;

import com.google.common.base.Function;
import com.google.common.collect.ImmutableList;
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
	protected Iterator<DirectedEvidence> perChrIterator(String chr) {
		FileSystemContext fsc = processContext.getFileSystemContext();
		return (Iterator<DirectedEvidence>)(Iterator)iterator( // TODO: is there a less dirty way to do this downcast?
				fsc.getBreakendVcfForChr(input, chr),
				fsc.getRealignmentBamForChr(input, chr));
	}
	@Override
	protected Iterator<DirectedEvidence> singleFileIterator() {
		FileSystemContext fsc = processContext.getFileSystemContext();
		return (Iterator<DirectedEvidence>)(Iterator)iterator(
			fsc.getBreakendVcf(input),
			fsc.getRealignmentBam(input));
	}
	private Iterator<VariantContextDirectedEvidence> iterator(File vcf, File realignment) {
		Iterator<SAMRecord> realignedIt = ImmutableList.<SAMRecord>of().iterator(); 
		if (isRealignmentComplete()) {
			realignedIt = new AutoClosingIterator<SAMRecord>(processContext.getSamReaderIterator(processContext.getSamReader(vcf)));
		}
		Iterator<VariantContextDirectedEvidence> it = new VariantContextDirectedEvidenceIterator(
				processContext,
				this,
				new AutoClosingIterator<VariantContext>(new VCFFileReader(vcf).iterator()),
				realignedIt);
		it = new DirectEvidenceWindowedSortingIterator<VariantContextDirectedEvidence>(processContext, 3 * this.getMetrics().getMaxFragmentSize(), it);
		return it;
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
		List<Closeable> toClose = Lists.newArrayList();
		SAMSequenceDictionary dict = processContext.getReference().getSequenceDictionary();
		try {
			if (processContext.shouldProcessPerChromosome()) {
				for (int i = 0; i < dict.size(); i++) {
					String seq = dict.getSequence(i).getSequenceName();
					List<Iterator<DirectedEvidence>> toMerge = Lists.newArrayList();
					for (SAMEvidenceSource bam : source) {
						Iterator<DirectedEvidence> it = bam.iterator(seq);
						if (toClose instanceof Closeable) toClose.add((Closeable)it);
						toMerge.add(it);
					}
					Iterator<DirectedEvidence> merged = Iterators.mergeSorted(toMerge, DirectedEvidenceOrder.ByNatural);
					assemble(merged, processContext.getFileSystemContext().getBreakendVcfForChr(input, seq), processContext.getFileSystemContext().getRealignmentFastqForChr(input, seq));
					for (Iterator<DirectedEvidence> x : toMerge) {
						CloserUtil.close(x);
					}
				}
			} else {
				List<Iterator<DirectedEvidence>> toMerge = Lists.newArrayList();
				for (SAMEvidenceSource bam : source) {
					Iterator<DirectedEvidence> it = bam.iterator();
					if (toClose instanceof Closeable) toClose.add((Closeable)it);
					toMerge.add(it);
				}
				Iterator<DirectedEvidence> merged = Iterators.mergeSorted(toMerge, DirectedEvidenceOrder.ByNatural);
				assemble(merged, processContext.getFileSystemContext().getBreakendVcf(input), processContext.getFileSystemContext().getRealignmentFastq(input));
				for (Iterator<DirectedEvidence> x : toMerge) {
					CloserUtil.close(x);
				}
			}
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
	private void processAssembly(Iterable<VariantContextDirectedEvidence> evidenceList, FastqBreakpointWriter fastqWriter, VariantContextWriter vcfWriter) {
    	if (evidenceList != null) {
	    	for (VariantContextDirectedEvidence a : evidenceList) {
	    		processAssembly(a, fastqWriter, vcfWriter);
	    	}
    	}
    }
	private void processAssembly(VariantContextDirectedEvidence evidence, FastqBreakpointWriter fastqWriter, VariantContextWriter vcfWriter) {
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
