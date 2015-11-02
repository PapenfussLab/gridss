package au.edu.wehi.idsv.alignment;

import java.io.File;
import java.io.IOException;
import java.lang.ProcessBuilder.Redirect;
import java.util.List;
import java.util.stream.Collectors;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Log;

public class ExternalProcessFastqAligner implements FastqAligner {
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
		log.info("Invoking external aligner");
		log.info(commandline.stream().collect(Collectors.joining(" ")));
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
		if (aligner.isAlive()) {
			log.error("External process still alive after alignment");
		}
		if (aligner.exitValue() != 0) {
			String msg = String.format("Subprocess terminated with with exit status %1$d. Alignment failed for %2$s", aligner.exitValue(), fastq);
			log.error(msg);
			throw new RuntimeException(msg);
		}
	}
}
