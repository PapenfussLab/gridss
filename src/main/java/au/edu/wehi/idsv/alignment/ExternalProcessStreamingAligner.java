package au.edu.wehi.idsv.alignment;

import java.io.BufferedOutputStream;
import java.io.Closeable;
import java.io.File;
import java.io.Flushable;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.ProcessBuilder.Redirect;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.Log;

/**
 * Performs alignment of the given records using an external alignment tools.
 * stdin and stdout of the alignment tool are hooked up 
 * 
 * @author Daniel Cameron
 *
 */
public class ExternalProcessStreamingAligner implements Closeable, Flushable {
	private static final int POLL_INTERVAL = 100;
	private static final Log log = Log.getInstance(ExternalProcessStreamingAligner.class);	
	private final AtomicInteger outstandingReads = new AtomicInteger(0);
	private final BlockingQueue<SAMRecord> buffer = new LinkedBlockingQueue<>();
	private final List<String> args;
	private final SamReaderFactory readerFactory;
	private Process aligner = null;
	private BasicFastqWriter toExternalProgram = null;
	private Thread reader = null;
	// The following are only needed for pretty error messages
	private final String commandlinestr;
	private final File reference;
	public ExternalProcessStreamingAligner(final SamReaderFactory readerFactory, final List<String> commandline, final File reference, final int threads) {
		this.readerFactory = readerFactory;
		this.reference = reference;
		this.args = commandline.stream()
				.map(s -> String.format(s, "-", reference.getAbsolutePath(), threads))
				.collect(Collectors.toList());
		this.commandlinestr = commandline.stream().collect(Collectors.joining(" "));
	}
	public synchronized void asyncAlign(FastqRecord fq) throws IOException {
		ensureAligner();
		outstandingReads.incrementAndGet();
		toExternalProgram.write(fq);
		toExternalProgram.flush();
	}
	public void ensureAligner() throws IOException {
		if (aligner == null) {
			log.info("Starting external aligner");
			log.info(commandlinestr);
			aligner = new ProcessBuilder(args)
					.redirectInput(Redirect.PIPE)
					.redirectOutput(Redirect.PIPE)
					.redirectError(Redirect.INHERIT)
					.start();
			toExternalProgram = new BasicFastqWriter(new PrintStream(new BufferedOutputStream(aligner.getOutputStream())));
			reader = new Thread(() -> readAllAlignments(readerFactory));
			reader.start();
		}
	}
	/***
	 * Flushes outstanding alignment requests.
	 * 
	 * Note that both external programs and the OS buffer both the input and output streams.
	 * To guarantee that the reads have been flushed, this methods closes the external
	 * program and as such, is a very expensive operation.
	 * @throws IOException 
	 * 
	 */
	public void flush() throws IOException {
		close();
	}
	/**
	 * Returns true if there is at least one outstanding alignment completed   
	 * @return
	 */
	public boolean hasAlignmentRecord() {
		return buffer.size() > 0;
	}
	/**
	 * Gets an alignment completed by the caller.
	 * 
	 * @return
	 * @throws IllegalStateException thrown when no alignment record is available from the aligner. Check if a record is available using hasAlignmentRecord() 
	 */
	public SAMRecord getAlignment() {
		if (!hasAlignmentRecord()) {
			throw new IllegalStateException("No alignments available. getAlignment() should only be called if at least one alignment record is available.");
		}
		SAMRecord r = buffer.poll();
		if (r == null) {
			throw new IllegalStateException("Empty buffer after flushing.");
		}
		return r;
	}
	private void readAllAlignments(final SamReaderFactory readerFactory) {
		SamReader fromExternalProgram = readerFactory.open(SamInputResource.of(aligner.getInputStream()));
		SAMRecordIterator it = fromExternalProgram.iterator();
		while (it.hasNext()) {
			SAMRecord r = it.next();
			outstandingReads.decrementAndGet();
			buffer.add(r);
		}
	}
	public void close() throws IOException {
		if (aligner == null) {
			// nothing to do
			return;
		}
		log.info("Waiting for external aligner to complete all alignments.");
		toExternalProgram.close();
		// and just to be sure we don't hit any more htsjdk bugs where they don't close the underlying stream
		aligner.getOutputStream().close(); 
		// wait for the aligner to complete all outstanding alignments
		// This doesn't deadlock as buffer is unbounded in size so we're guaranteed to be able to
		// read the entire output stream without blocking
		while (outstandingReads.get() > 0) {
			try {
				Thread.sleep(POLL_INTERVAL);
			} catch (InterruptedException e) {
				log.warn(e);
				return;
			}
		}
		// reader thread will have completed when it hits then end of the output stream 
		ExternalProcessHelper.shutdownAligner(aligner, commandlinestr, reference);
		log.info("External alignments complete");
		aligner = null;
		reader = null;
		toExternalProgram = null;
	}
}
