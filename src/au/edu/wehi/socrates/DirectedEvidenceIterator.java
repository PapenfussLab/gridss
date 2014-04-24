package au.edu.wehi.socrates;

import java.io.File;
import java.util.ArrayDeque;
import java.util.Collection;
import java.util.PriorityQueue;
import java.util.Queue;

import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.fastq.FastqWriterFactory;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;

import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;

import au.edu.wehi.socrates.util.SAMRecordSummary;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

/**
 * Iterates through evidence for a breakpoint at 
 * @author Daniel Cameron
 *
 */
public class DirectedEvidenceIterator extends AbstractIterator<DirectedEvidence> {
	private static final int INITIAL_BUFFER_SIZE = 1024;
	private final PeekingIterator<SAMRecord> svIterator;
	private final PeekingIterator<VariantContext> vcfIterator;
	private final LinearGenomicCoordinate linear;
	private final SequentialNonReferenceReadPairFactory mateFactory;
	private final SequentialRealignedBreakpointFactory realignFactory;
	private final PriorityQueue<DirectedEvidence> calls = new PriorityQueue<DirectedEvidence>(INITIAL_BUFFER_SIZE, new DirectedEvidenceStartCoordinateComparator());
	private final int maxFragmentSize;
	private long currentLinearPosition = -1;
	/**
	 * Iterates over breakpoint evidence
	 * @param sv structural variation supporting reads. Must be sorted by coordinate
	 * @param mate mate reads of structural variation supporting reads. Must be sorted by mate coordinate
	 * @param realign breakpoint realignments. Must be sorted by source @see DirectedBreakpoint start coordinate 
	 * @param vcf Assembly-based @see DirectedBreakpoint. Must be sorted by coordinate
	 */
	public DirectedEvidenceIterator(
			PeekingIterator<SAMRecord> sv,
			PeekingIterator<SAMRecord> mate,
			PeekingIterator<SAMRecord> realign,
			PeekingIterator<VariantContext> vcf,
			SAMSequenceDictionary dictionary,
			int maxFragmentSize) {
		this.linear = new LinearGenomicCoordinate(dictionary);
		this.svIterator = sv;
		this.vcfIterator = vcf;
		this.mateFactory = new SequentialNonReferenceReadPairFactory(mate);
		this.realignFactory = new SequentialRealignedBreakpointFactory(realign);
		this.maxFragmentSize = maxFragmentSize;
	}
	@Override
	protected DirectedEvidence computeNext() {
		do {
			if (!calls.isEmpty() && linear.getStartLinearCoordinate(calls.peek().getBreakpointLocation()) < currentLinearPosition - maxFragmentSize) {
				DirectedEvidence evidence = calls.poll();
				// Match realigned at output time since input BAM processing is in order of alignment start,
				// not putative breakpoint position 
				if (evidence.getClass() == DirectedBreakpoint.class) {
					DirectedBreakpoint bp = (DirectedBreakpoint)evidence;
					bp.setRealigned(realignFactory.findRealignedSAMRecord(bp));
				}
				return evidence;
			}
		} while (advance());
		return endOfData();
	}
	/**
	 * Advances the input streams to the next
	 * @return true if records were consumed, false if no more input exists
	 */
	private boolean advance() {
		long nextPos = Long.MAX_VALUE;
		if (svIterator != null && svIterator.hasNext()) {
			nextPos = Math.min(nextPos, linear.getLinearCoordinate(svIterator.peek().getReferenceIndex(), svIterator.peek().getAlignmentStart()));
		}
		if (vcfIterator != null && vcfIterator.hasNext()) {
			nextPos = Math.min(nextPos, linear.getLinearCoordinate(vcfIterator.peek().getChr(), vcfIterator.peek().getStart()));
		}
		// no records
		if (nextPos == Long.MAX_VALUE) return false;
		while (svIterator != null && svIterator.hasNext() && linear.getLinearCoordinate(svIterator.peek().getReferenceIndex(), svIterator.peek().getAlignmentStart()) == nextPos) {
			processRead(svIterator.next());
		}
		while (vcfIterator != null && vcfIterator.hasNext() && linear.getLinearCoordinate(vcfIterator.peek().getChr(), vcfIterator.peek().getStart()) == nextPos) {
			processVariant(vcfIterator.next());
		}
		currentLinearPosition = nextPos;
		return true;
	}
	private void processVariant(VariantContext variant) {
		AssemblyEvidence evidence = new AssemblyEvidence(variant);
		if (evidence.isValid()) {
			calls.add(evidence);
		}
	}
	private void processRead(SAMRecord record) {
		if (SAMRecordSummary.getStartSoftClipLength(record) > 0) {
			calls.add(new SoftClipEvidence(BreakpointDirection.Backward, record));
		}
		if (SAMRecordSummary.getEndSoftClipLength(record) > 0) {
			calls.add(new SoftClipEvidence(BreakpointDirection.Forward, record));
		}
		if (SAMRecordSummary.isPartOfNonReferenceReadPair(record)) {
			NonReferenceReadPair pair = mateFactory.createNonReferenceReadPair(record, maxFragmentSize);
			if (pair != null) {
				calls.add(pair);
			}
		}
	}
}
