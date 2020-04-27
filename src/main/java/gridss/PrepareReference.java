package gridss;

import au.edu.wehi.idsv.GenomicProcessingContext;
import au.edu.wehi.idsv.alignment.BwaAligner;
import au.edu.wehi.idsv.picard.TwoBitBufferedReferenceSequenceFile;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;

import java.io.File;
import java.util.Locale;

@CommandLineProgramProperties(
        summary = "Creates index and cache files required by GRIDSS.",
        oneLineSummary = "Creates index and cache files required by GRIDSS.",
        programGroup = gridss.cmdline.programgroups.DataConversion.class
)
public class PrepareReference extends CommandLineProgram {
	private static final Log log = Log.getInstance(PrepareReference.class);
	@Argument(doc="Creates sequence dictionary.", optional=true)
	public boolean CREATE_SEQUENCE_DICTIONARY = true;
	@Argument(doc="Creates 2-bit encoded reference cache file used by GRIDSS.", optional=true)
	public boolean CREATE_GRIDSS_REFERENCE_CACHE = true;
	@Argument(doc="Creates bwa index file used by the in-process bwa aligner.", optional=true)
	public boolean CREATE_BWA_INDEX_IMAGE = true;
    @Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
    	Locale.setDefault(Locale.ROOT);
		try {
			validateParameters();
			if (CREATE_SEQUENCE_DICTIONARY) {
				if (!ReferenceCommandLineProgram.ensureSequenceDictionary(REFERENCE_SEQUENCE)) {
					log.info("Sequence dictionary found.");
				}
			}
			File cache = GenomicProcessingContext.getGridssCacheFileForReference(REFERENCE_SEQUENCE);
			if (CREATE_GRIDSS_REFERENCE_CACHE) {
				if (!cache.exists()) {
					log.info("Creating GRIDSS reference cache file " + cache);
					ReferenceSequenceFile ref = new IndexedFastaSequenceFile(REFERENCE_SEQUENCE);
					TwoBitBufferedReferenceSequenceFile tbbrsf = new TwoBitBufferedReferenceSequenceFile(ref);
					tbbrsf.save(cache);
				} else {
					log.info("Found " + cache);
				}
			}
			File bwaImage = BwaAligner.getBwaIndexFileFor(REFERENCE_SEQUENCE);
			if (CREATE_BWA_INDEX_IMAGE) {
				if (!bwaImage.exists() || bwaImage.length() == 0) {
					if (bwaImage.exists()) {
						log.info("Removing empty file " + cache);
						bwaImage.delete();
					}
					log.info("Creating BWA index image file " + cache);
					BwaAligner.createBwaIndexFor(REFERENCE_SEQUENCE);
				} else {
					log.info("Found " + bwaImage);
				}
			}
		} catch (Exception e) {
			log.error(e);
			return -1;
		}
    	return 0;
	}
    
	private void validateParameters() {
    	IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
	}

	@Override
	protected String[] customCommandLineValidation() {
    	if (REFERENCE_SEQUENCE == null) {
    		return new String[] { "REFERENCE_SEQUENCE is required."};
		}
		if (!REFERENCE_SEQUENCE.isFile()) {
			return new String[] { "REFERENCE_SEQUENCE file not found."};
		}
		return super.customCommandLineValidation();
	}

	public static void main(String[] argv) {
        System.exit(new PrepareReference().instanceMain(argv));
    }
}
