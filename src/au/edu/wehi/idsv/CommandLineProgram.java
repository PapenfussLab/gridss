package au.edu.wehi.idsv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;

import java.io.File;
import java.io.IOException;

import au.edu.wehi.idsv.vcf.VcfConstants;
import picard.cmdline.Option;

public abstract class CommandLineProgram extends picard.cmdline.CommandLineProgram {
	@Option(doc = "Breakends are written to VCF files as VCF v4.1 compatible breakpoints to a placeholder contig " + VcfConstants.VCF41BREAKEND_REPLACEMENT,
            optional = true,
            shortName = "VCF41")
	public boolean VCF41_COMPATIBLE = true;
	public SamReaderFactory getSamReaderFactory() {
		return SamReaderFactory.makeDefault()
				.validationStringency(ValidationStringency.LENIENT);
	}
	public File ensureFileExists(final File file) {
    	if (!file.exists()) {
    		throw new RuntimeException(String.format("Required input file %s is missing. Have the earlier pipeline commands executed successfully?", file));
    	}
    	return file;
    }
	public ProcessingContext getContext(File reference, File input) throws IOException {
		final ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(reference);
		final SAMSequenceDictionary dictionary = ref.getSequenceDictionary();
		if (dictionary == null) {
			throw new IllegalArgumentException(String.format("Missing .dict for %s. Creating using picard tools CreateSequenceDictionary.", reference));
		}
    	final RelevantMetrics metrics = new RelevantMetrics(FileNamingConvention.getMetrics(input));
    	final ProcessingContext processContext = new ProcessingContext(ref, dictionary, metrics);
    	processContext.setVcf41Mode(VCF41_COMPATIBLE);
    	return processContext;
	}
}
