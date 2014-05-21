package au.edu.wehi.socrates;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;

import java.io.File;
import java.io.IOException;

public abstract class CommandLineProgram extends picard.cmdline.CommandLineProgram {
	public SamReaderFactory getSamReaderFactory() {
		return SamReaderFactory.makeDefault()
				.validationStringency(ValidationStringency.LENIENT);
	}
	public static File ensureFileExists(final File file) {
    	if (!file.exists()) {
    		throw new RuntimeException(String.format("Required input file %s is missing. Have the earlier pipeline commands executed successfully?", file));
    	}
    	return file;
    }
	public static ProcessingContext getContext(File reference, File input) throws IOException {
		final ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(reference);
		final SAMSequenceDictionary dictionary = ref.getSequenceDictionary();
		if (dictionary == null) {
			throw new IllegalArgumentException(String.format("Missing .dict for ", reference, ". Creating using picard tools CreateSequenceDictionary."));
		}
    	final RelevantMetrics metrics = new RelevantMetrics(FileNamingConvention.getMetrics(input));
    	final ProcessingContext processContext = new ProcessingContext(ref, dictionary, metrics);
    	return processContext;
	}
}
