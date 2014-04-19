package au.edu.wehi.socrates;

import java.util.List;

import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;

import com.google.common.base.Charsets;

import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMUtils;

public class SAMRecordEvidenceProcessorOrchestrator {
	public static void process(
			SequentialSAMRecordBuffer itbuffer,
			List<SequentialDirectedBreakpointEvidenceProcessor> processors,
			int maxFragmentSize,
			FastqWriter fastqWriter,
			VariantContextWriter vcfWriter) {
		int position = 1;
		while (itbuffer.hasReads()) {
			// using VCF breakpoint position definition
			// Forward evidence
			for (SAMRecord r = itbuffer.getNextEndingSoftClippedRead(position); r != null; r = itbuffer.getNextEndingSoftClippedRead(position)) {
				// add to forward working set
				for (SequentialDirectedBreakpointEvidenceProcessor processor : processors) {
					processor.addSoftclipSupportingForwardBreakpoint(r);
				}
			}
			for (NonReferenceReadPair r = itbuffer.getNextNonReferenceForwardReadPair(position); r != null; r = itbuffer.getNextNonReferenceForwardReadPair(position)) {
				// add to forward working set
				for (SequentialDirectedBreakpointEvidenceProcessor processor : processors) {
					processor.addReadPairSupportingForwardBreakpoint(r);
				}
			}
			// Backward evidence
			for (SAMRecord r = itbuffer.getNextStartingSoftClippedRead(position + 1); r != null; r = itbuffer.getNextStartingSoftClippedRead(position + 1)) {
				// add to backward working set
				for (SequentialDirectedBreakpointEvidenceProcessor processor : processors) {
					processor.addSoftclipSupportingBackwardBreakpoint(r);
				}
			}
			for (NonReferenceReadPair r = itbuffer.getNextNonReferenceBackwardReadPair(position + maxFragmentSize); r != null; r = itbuffer.getNextNonReferenceBackwardReadPair(position + maxFragmentSize)) {
				// add to backward working set
				for (SequentialDirectedBreakpointEvidenceProcessor processor : processors) {
					processor.addReadPairSupportingBackwardBreakpoint(r);
				}
			}
			// process genomic position
			for (SequentialDirectedBreakpointEvidenceProcessor processor : processors) {
				for (DirectedBreakpoint bp : processor.getEvidenceAtPosition(position)) {
					byte[] sequence = bp.getBreakpointSequence();
					if (sequence != null) {
						FastqRecord fq = new FastqRecord(
								bp.getBreakpointID(),
								new String(sequence, Charsets.US_ASCII),
								"",
								SAMUtils.phredToFastq(bp.getBreakpointQuality()));
						fastqWriter.write(fq);
					}
					VariantContext variant = bp.getVariantContext();
					if (variant != null) {
						vcfWriter.add(variant);
					}
				}
			}
			// advance genomic position to next callable position
			position++;
		}
	}
}
