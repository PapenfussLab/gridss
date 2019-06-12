package au.edu.wehi.idsv.alignment;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;

import java.io.IOException;

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
	 * Number of fastq records that have been aligned but not yet returned by getAlignment()
	 * @return
	 */
	int processedAlignmentRecords();
	
	/**
	 * Number of fastq records sent to the aligner for which no response has yet been received
	 * @return
	 */
	int outstandingAlignmentRecord();

	/**
	 * Gets an alignment completed by the caller.
	 * 
	 * @return
	 * @throws IllegalStateException thrown when no alignment record is available from the aligner. Check if a record is available using hasAlignmentRecord() 
	 */
	SAMRecord getAlignment();

	void close() throws IOException;
}