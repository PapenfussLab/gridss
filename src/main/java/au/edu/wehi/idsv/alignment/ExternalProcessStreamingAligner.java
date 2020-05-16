package au.edu.wehi.idsv.alignment;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.MessageThrottler;
import htsjdk.samtools.*;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.Log;
import org.apache.commons.lang3.SystemUtils;

import java.io.*;
import java.lang.ProcessBuilder.Redirect;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

/**
 * Performs alignment of the given records using an external alignment tools.
 * stdin and stdout of the alignment tool are hooked up.
 * 
 * Iterator methods block until an input record has been aligned. The StreamingAligner
 * interface allows for more fine-grain control over record processing.
 * 
 * @author Daniel Cameron
 *
 */
public class ExternalProcessStreamingAligner implements Closeable, Flushable, StreamingAligner {
	private static final int POLL_INTERVAL = 1000;
	private static final int OUTPUT_BUFFER_SIZE = 1024;
	private static final Log log = Log.getInstance(ExternalProcessStreamingAligner.class);
	private final AtomicInteger outstandingReads = new AtomicInteger(0);
	private final BlockingQueue<SAMRecord> buffer = new ArrayBlockingQueue<>(OUTPUT_BUFFER_SIZE);
	private final List<String> args;
	private final SamReaderFactory readerFactory;
	private final SAMSequenceDictionary dict;
	private Process aligner = null;
	private BasicFastqWriter toExternalProgram = null;
	private Thread reader = null;
	// The following are only needed for pretty error messages
	private final String commandlinestr;
	private final File reference;
	public ExternalProcessStreamingAligner(final SamReaderFactory readerFactory, final List<String> commandline, final File reference, final int threads, final SAMSequenceDictionary dict) {
		this.readerFactory = readerFactory;
		this.reference = reference;
		this.dict = dict;
		this.args = commandline.stream()
				.map(s -> String.format(s, "-", reference.getPath(), threads))
				.collect(Collectors.toList());
		this.commandlinestr = args.stream().collect(Collectors.joining(" "));
	}
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
			List<String> commandline = args;
			if (SystemUtils.IS_OS_WINDOWS) {
				// WSL path conversion
				commandline = commandline.stream()
						.map(s -> s.replace('\\', '/'))
						.collect(Collectors.toList());
			}
			aligner = new ProcessBuilder(commandline)
					.redirectInput(Redirect.PIPE)
					.redirectOutput(Redirect.PIPE)
					.redirectError(Redirect.INHERIT)
					.start();
			toExternalProgram = new BasicFastqWriter(new PrintStream(new BufferedOutputStream(aligner.getOutputStream())));
			reader = new Thread(() -> readAllAlignments(readerFactory));
			reader.setName("ExternalProcessStreamingAligner");
			reader.start();
		}
	}
	@Override
	public void flush() throws IOException {
		if (aligner != null) {
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
				}
			}
			// reader thread will have completed when it hits then end of the output stream
			ExternalProcessHelper.shutdownAligner(aligner, commandlinestr, reference, null);
			log.info("External alignments complete");
			try {
				reader.join();
			} catch (InterruptedException e) {
				log.warn(e);
			}
		}
		aligner = null;
		reader = null;
		toExternalProgram = null;
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
		SAMRecord r = buffer.poll();
		if (r == null) {
			throw new IllegalStateException("No alignments available. getAlignment() should only be called if at least one alignment record is available.");
		}
		return r;
	}
	private void readAllAlignments(final SamReaderFactory readerFactory) {
		try {
			SamReader fromExternalProgram = readerFactory.open(SamInputResource.of(aligner.getInputStream()));
			SAMRecordIterator it = fromExternalProgram.iterator();
			while (it.hasNext()) {
				SAMRecord r = it.next();
				if (SAMRecordUtil.forceValidContigBounds(r, dict)) {
					if (!MessageThrottler.Current.shouldSupress(log, "streaming aligner out of bounds")) {
						log.warn(String.format("Streamed aligner returned out of bounds alignment. %s adjusted to %s:%d %s", dict.getSequence(r.getReferenceIndex()).getSequenceName(), r.getAlignmentStart(), r.getCigarString()));
					}
				}
				buffer.put(r);
				outstandingReads.decrementAndGet();
			}
			log.info(String.format("Reader thread complete. %s reads in output buffer", buffer.size()));
		} catch (InterruptedException ie) {
			log.warn(ie, "reader thread interrupted");
		}
	}

	@Override
	public void close() throws IOException {
		flush();
	}
}
