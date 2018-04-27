package au.edu.wehi.idsv.alignment;

import java.io.BufferedOutputStream;
import java.io.Closeable;
import java.io.File;
import java.io.Flushable;
import java.io.IOException;
import java.lang.ProcessBuilder.Redirect;
import java.util.Iterator;
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
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.NonFlushingBasicFastqWriter;
import htsjdk.samtools.util.Log;

/**
 * Performs alignment of the given records using an external alignment tools.
 * stdin and stdout of the alignment tool are hooked up 
 * 
 * @author Daniel Cameron
 *
 */
public class ExternalProcessStreamingAligner implements Closeable, Flushable, StreamingAligner, Iterator<SAMRecord> {
	private static final int POLL_INTERVAL = 1000;
	private static final Log log = Log.getInstance(ExternalProcessStreamingAligner.class);	
	private final AtomicInteger outstandingReads = new AtomicInteger(0);
	private final BlockingQueue<SAMRecord> buffer = new LinkedBlockingQueue<>();
	private final List<String> args;
	private final SamReaderFactory readerFactory;
	private Process aligner = null;
	private NonFlushingBasicFastqWriter toExternalProgram = null;
	private Thread reader = null;
	// The following are only needed for pretty error messages
	private final String commandlinestr;
	private final File reference;
	public ExternalProcessStreamingAligner(final SamReaderFactory readerFactory, final List<String> commandline, final File reference, final int threads) {
		this.readerFactory = readerFactory;
		this.reference = reference;
		this.args = commandline.stream()
				.map(s -> String.format(s, "-", reference.getPath(), threads))
				.collect(Collectors.toList());
		this.commandlinestr = args.stream().collect(Collectors.joining(" "));
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.alignment.StreamingAligner#asyncAlign(htsjdk.samtools.fastq.FastqRecord)
	 */
	@Override
	public synchronized void asyncAlign(FastqRecord fq) throws IOException {
		ensureAligner();
		outstandingReads.incrementAndGet();
		toExternalProgram.write(fq);
		toExternalProgram.flush();
	}
	private void ensureAligner() throws IOException {
		if (aligner == null) {
			log.info("Starting external aligner");
			log.info(commandlinestr);
			aligner = new ProcessBuilder(args)
					.redirectInput(Redirect.PIPE)
					.redirectOutput(Redirect.PIPE)
					.redirectError(Redirect.INHERIT)
					.start();
			toExternalProgram = new NonFlushingBasicFastqWriter(new BufferedOutputStream(aligner.getOutputStream()));
			reader = new Thread(() -> readAllAlignments(readerFactory));
			reader.setName("ExternalProcessStreamingAligner");
			reader.start();
		}
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.alignment.StreamingAligner#flush()
	 */
	@Override
	public void flush() throws IOException {
		close();
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.alignment.StreamingAligner#hasAlignmentRecord()
	 */
	@Override
	public boolean hasAlignmentRecord() {
		return buffer.size() > 0;
	}
	@Override
	public int processedAlignmentRecords() {
		return buffer.size();
	}
	@Override
	public int outstandingAlignmentRecord() {
		return outstandingReads.get();
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.alignment.StreamingAligner#getAlignment()
	 */
	@Override
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
			buffer.add(r);
			outstandingReads.decrementAndGet();
		}
		log.info(String.format("Reader thread complete. %s reads in output buffer", buffer.size()));
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.alignment.StreamingAligner#close()
	 */
	@Override
	public synchronized void close() throws IOException {
		if (aligner == null) {
			// nothing to do
			return;
		}
		log.info("Waiting for external aligner to complete all alignments.");
		toExternalProgram.flush();
		aligner.getOutputStream().flush();
		toExternalProgram.close();
		// and just to be sure we don't hit any more htsjdk bugs where they don't close the underlying stream
		aligner.getOutputStream().close();
		// wait for the aligner to complete all outstanding alignments
		// This doesn't deadlock as buffer is unbounded in size so we're guaranteed to be able to
		// read the entire output stream without blocking
		while (outstandingReads.get() > 0) {
			try {
				log.debug(String.format("%d alignments outstanding", outstandingReads.get()));
				Thread.sleep(POLL_INTERVAL);
			} catch (InterruptedException e) {
				log.warn(e);
				return;
			}
		}
		// reader thread will have completed when it hits then end of the output stream 
		ExternalProcessHelper.shutdownAligner(aligner, commandlinestr, reference);
		log.info("External alignments complete");
		try {
			reader.join();
		} catch (InterruptedException e) {
			log.warn(e);
		}
		aligner = null;
		reader = null;
		toExternalProgram = null;
	}
	@Override
	public boolean hasNext() {
		return hasAlignmentRecord();
	}
	@Override
	public SAMRecord next() {
		return getAlignment();
	}
}
