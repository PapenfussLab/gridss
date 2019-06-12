package au.edu.wehi.idsv.alignment;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Log;
import org.apache.commons.lang3.SystemUtils;

import java.io.File;
import java.io.IOException;
import java.lang.ProcessBuilder.Redirect;
import java.util.List;
import java.util.stream.Collectors;

public class ExternalProcessFastqAligner implements FastqAligner {
	/**
	 * Number of seconds to wait for the external aligner to shut down
	 */
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
				.map(s -> String.format(s, fastq.getPath(), reference.getPath(), threads))
				.collect(Collectors.toList());
		if (SystemUtils.IS_OS_WINDOWS) {
			// WSL path conversion
			commandline = commandline.stream()
				.map(s -> s.replace('\\', '/'))
				.collect(Collectors.toList());
		}
		String commandlinestr = commandline.stream().collect(Collectors.joining(" "));
		log.info("Invoking external aligner");
		log.info(commandlinestr);
		Process aligner = new ProcessBuilder(commandline)
				.redirectError(Redirect.INHERIT)
				.start();
		final SamReader reader = readerFactory.open(SamInputResource.of(aligner.getInputStream()));
		final SAMFileHeader header = reader.getFileHeader();
		try (final SAMFileWriter writer = writerFactory.clone().setCompressionLevel(0).makeWriter(header, false, output, reference)) {
			final SAMRecordIterator it = reader.iterator();
			while (it.hasNext()) {
				writer.addAlignment(it.next());
			}
		}
		ExternalProcessHelper.shutdownAligner(aligner, commandlinestr, reference);
	}
}
