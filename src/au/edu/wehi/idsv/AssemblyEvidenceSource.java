package au.edu.wehi.idsv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.util.Iterator;
import java.util.List;

import au.edu.wehi.idsv.debruijn.anchored.DeBruijnAnchoredAssembler;
import au.edu.wehi.idsv.debruijn.subgraph.DeBruijnSubgraphAssembler;

import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;


public class AssemblyEvidenceSource extends EvidenceSource {
	private static final Log log = Log.getInstance(SAMEvidenceSource.class);
	private final List<SAMEvidenceSource> source;
	/**
	 * Generates assembly evidence based on the given evidence
	 * @param evidence evidence for creating assembly
	 * @param intermediateFileLocation location to store intermediate files
	 */
	public AssemblyEvidenceSource(ProcessingContext processContext, List<SAMEvidenceSource> evidence, File intermediateFileLocation) {
		super(processContext, intermediateFileLocation);
		this.source = evidence;
	}
	@Override
	public Iterator<DirectedEvidence> iterator() {
		ensureAssembly();
		if (!assemblyDone()) {
			throw new IllegalStateException(String.format("Missing intermediate files for ", input));
		} else if (!isRealignmentComplete()) {
			log.debug("Missing realignment bam. Traversing breakends only.");
		}
		return null;
	}
	@Override
	public Iterator<DirectedEvidence> iterator(String chr) {
		// TODO Auto-generated method stub
		return null;
	}
	private void ensureAssembly() {
		if (!assemblyDone()) {
			log.info("START evidence assembly ", input);
			assemble();
			log.info("SUCCESS evidence assembly ", input);
		}
	}
	private boolean assemblyDone() {
		boolean done = true;
		FileSystemContext fsc = processContext.getFileSystemContext();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				done |= checkIntermediate(fsc.getBreakendVcfForChr(input, seq.getSequenceName()));
				done |= checkIntermediate(fsc.getRealignmentFastqForChr(input, seq.getSequenceName()));
			}
		} else {
			done |= checkIntermediate(fsc.getBreakendVcf(input));
			done |= checkIntermediate(fsc.getRealignmentFastq(input));
		}
		return done;
	}
	private void assemble() {
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
					Iterator<DirectedEvidence> merged = Iterators.mergeSorted(toMerge, DirectedEvidenceOrder.ByStartEnd);
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
				Iterator<DirectedEvidence> merged = Iterators.mergeSorted(toMerge, DirectedEvidenceOrder.ByStartEnd);
				assemble(merged, processContext.getFileSystemContext().getBreakendVcf(input));
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
	private void assemble(Iterator<DirectedEvidence> it, File breakendVcf) {
		ReadEvidenceAssembler assembler = getAssembler();
		while (it.hasNext()) {
			DirectedEvidence readEvidence = it.next();
			progress.record(readEvidence.)
			// Need to process assembly evidence first since assembly calls are made when the
			// evidence falls out of scope so processing a given position will emit evidence
			// for a previous position (for it is only at this point we know there is no more
			// evidence for the previous position).
			processAssemblyEvidence(assembler.addEvidence(readEvidence), fastqWriter, vcfWriter);
			/*if (readEvidence instanceof SoftClipEvidence) {
				SoftClipEvidence sce = (SoftClipEvidence)readEvidence;
				if (sce.getMappingQuality() >= MIN_MAPQ &&
						sce.getSoftClipLength() >= MIN_BREAKEND_REALIGN_LENGTH &&
						sce.getAlignedPercentIdentity() >= MIN_PERCENT_IDENTITY &&
						sce.getAverageClipQuality() >= MIN_LONG_SC_BASE_QUALITY) {
					fastqWriter.write(sce);
				}
			}*/
		}
		processAssemblyEvidence(assembler.endOfEvidence(), fastqWriter, vcfWriter);
	}
	private ReadEvidenceAssembler getAssembler() {
    	switch (processContext.getAssemblyParameters().method) {
	    	case DEBRUIJN_PER_POSITION:
	    		return new DeBruijnAnchoredAssembler(processContext);
	    	case DEBRUIJN_SUBGRAPH:
	    		return new DeBruijnSubgraphAssembler(processContext);
	    	default:
	    		throw new IllegalArgumentException("Unknown assembly method.");
    	}
    }
}
