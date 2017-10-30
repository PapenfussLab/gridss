package au.edu.wehi.idsv.alignment;

import java.io.IOException;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;

public interface StreamingAligner {

	void asyncAlign(FastqRecord fq) throws IOException;

	/***
	 * Flushes outstanding alignment requests.
	 * 
	 * Note that both external programs and the OS buffer both the input and output streams.
	 * To guarantee that the reads have been flushed, this methods closes the external
	 * program and as such, is a very expensive operation.
	 * @throws IOException 
	 * 
	 */
	void flush() throws IOException;

	/**
	 * Returns true if there is at least one outstanding alignment completed   
	 * @return
	 */
	boolean hasAlignmentRecord();

	/**
	 * Gets an alignment completed by the caller.
	 * 
	 * @return
	 * @throws IllegalStateException thrown when no alignment record is available from the aligner. Check if a record is available using hasAlignmentRecord() 
	 */
	SAMRecord getAlignment();

	void close() throws IOException;

}