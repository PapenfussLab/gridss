package au.edu.wehi.idsv.alignment;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.lang.ProcessBuilder.Redirect;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

public class ExternalProcessFastqAligner implements FastqAligner {
	/**
	 * Number of seconds to wait for the external aligner to shut down
	 */
	private static final int SHUTDOWN_GRACE_PERIOD = 10;
	private static final Log log = Log.getInstance(ExternalProcessFastqAligner.class);
	private final List<String> template;
	private final SamReaderFactory readerFactory;
	private final SAMFileWriterFactory writerFactory;
	public ExternalProcessFastqAligner(
			final SamReaderFactory readerFactory,
			final SAMFileWriterFactory writerFactory,
			final List<String> commandline) {
		this.template = commandline;
		this.readerFactory = readerFactory;
		this.writerFactory = writerFactory;
	}
	@Override
	public void align(final File fastq, final File output, final File reference, final int threads) throws IOException {
		List<String> commandline = template.stream()
				.map(s -> String.format(s, fastq.getAbsolutePath(), reference.getAbsolutePath(), threads))
				.collect(Collectors.toList());
		String commandlinestr = commandline.stream().collect(Collectors.joining(" "));
		log.info("Invoking external aligner");
		log.info(commandlinestr);
		Process aligner = new ProcessBuilder(commandline)
				.redirectError(Redirect.INHERIT)
				.directory(output.getParentFile())
				.start();
		final SamReader reader = readerFactory.open(SamInputResource.of(aligner.getInputStream()));
		final SAMFileHeader header = reader.getFileHeader();
		try (final SAMFileWriter writer = writerFactory.makeWriter(header, false, output, reference)) {
			final SAMRecordIterator it = reader.iterator();
			while (it.hasNext()) {
				writer.addAlignment(it.next());
			}
		}
		try {
			aligner.waitFor(SHUTDOWN_GRACE_PERIOD, TimeUnit.SECONDS);
		} catch (InterruptedException e) {
			// Restore the interrupted status
            Thread.currentThread().interrupt();
		}
		if (aligner.isAlive()) {
			log.error("External process still alive after alignment");
		}
		if (aligner.exitValue() != 0) {
			String msg = String.format("Subprocess terminated with with exit status %1$d. Alignment failed for %2$s executing \"%s\". Can you run the alignment command from the command line? Is the aligner on PATH? Did you build an index with prefix %3$s?", aligner.exitValue(), fastq, commandlinestr, reference);
			log.error(msg);
			throw new RuntimeException(msg);
		}
	}
}
