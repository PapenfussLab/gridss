package au.edu.wehi.socrates;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public abstract class CommandLineProgram extends picard.cmdline.CommandLineProgram {
	public SamReaderFactory getSamReaderFactory() {
		return SamReaderFactory.makeDefault()
				.validationStringency(ValidationStringency.LENIENT);
	}
}
