package au.edu.wehi.idsv;

import java.io.Closeable;
import java.util.Arrays;

import org.apache.commons.lang3.ArrayUtils;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.SequenceUtil;

/**
 * Writes breakpoint realignment records
 * @author Daniel Cameron
 *
 */
public class FastqBreakpointWriter implements Closeable {
	public FastqWriter backing;
	public FastqBreakpointWriter(FastqWriter writer) {
		this.backing = writer;
	}
	public void write(final DirectedEvidence breakpoint) {
		byte[] seq = breakpoint.getBreakendSequence();
		if (seq == null || seq.length == 0) throw new IllegalArgumentException(String.format("Breakend sequence unknown"));
		FastqRecord fastq = BreakpointFastqEncoding.getRealignmentFastq(breakpoint);
		backing.write(fastq);
	}
	/**
	 * Writes nested records for the given realigned breakend.
	 * If a realigned breakend contains soft clips, the breakend sequence may span multiple breakpoints.
	 * By realigning the breakend soft clips, we can find the alignment locations of all breakend bases.
	 * @param breakend breakend to realign
	 * @param minNestedLength
	 */
	public void writeRealignment(final SAMRecord breakend, int minRealignmentLength) {
		if (breakend.getReadUnmappedFlag()) return;
		int startClipLength = SAMRecordUtil.getStartSoftClipLength(breakend);
		int endClipLength = SAMRecordUtil.getEndSoftClipLength(breakend);
		if (startClipLength < minRealignmentLength && endClipLength < minRealignmentLength) return;
		
		String fastqid = breakend.getReadName();
		int refIndex = BreakpointFastqEncoding.getEncodedReferenceIndex(fastqid);
		int start = BreakpointFastqEncoding.getEncodedStartPosition(fastqid);
		int offset = BreakpointFastqEncoding.getEncodedBreakendOffset(fastqid);
		String id = BreakpointFastqEncoding.getEncodedID(fastqid);
		byte[] seq = breakend.getReadBases();
		byte[] qual = breakend.getBaseQualities();
		if (breakend.getReadNegativeStrandFlag()) {
			seq = Arrays.copyOf(seq, seq.length);
			SequenceUtil.reverseComplement(seq);
			qual = Arrays.copyOf(seq, seq.length);
			ArrayUtils.reverse(qual);
			int tmp = endClipLength;
			endClipLength = startClipLength;
			startClipLength = tmp;
		}
		if (startClipLength >= minRealignmentLength) {
			backing.write(BreakpointFastqEncoding.getRealignmentFastq(
					Arrays.copyOf(seq, startClipLength), 
					Arrays.copyOf(qual, startClipLength),
					refIndex,
					start,
					offset,
					id));
		}
		if (endClipLength >= minRealignmentLength) {
			backing.write(BreakpointFastqEncoding.getRealignmentFastq(
					Arrays.copyOfRange(seq, seq.length - endClipLength, seq.length), 
					Arrays.copyOfRange(qual, qual.length - endClipLength, qual.length),
					refIndex,
					start,
					offset + seq.length - endClipLength,
					id));
		}
	}
	public void close() {
		backing.close();
	}
}
