package au.edu.wehi.idsv;

import java.io.Closeable;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;

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
	public void close() {
		backing.close();
	}
}
