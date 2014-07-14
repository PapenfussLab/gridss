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
	public void write(final DirectedBreakpoint breakpoint) {
		FastqRecord fastq = BreakpointFastqEncoding.getRealignmentFastq(breakpoint);
		backing.write(fastq);
	}
	public void close() {
		backing.close();
	}
}
