package au.edu.wehi.socrates;

import java.io.File;
import java.io.IOException;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMSequenceDictionary;

import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;

/**
 * Writes breakpoint realignment records
 * @author Daniel Cameron
 *
 */
public class FastqBreakpointWriter {
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
