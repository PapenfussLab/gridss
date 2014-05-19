package au.edu.wehi.socrates;

import java.io.File;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

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
}
