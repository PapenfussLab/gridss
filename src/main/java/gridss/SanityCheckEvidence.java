package gridss;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.configuration.GridssConfiguration;
import au.edu.wehi.idsv.validation.PairedEvidenceTracker;
import gridss.cmdline.FullEvidenceCommandLineProgram;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import java.io.IOException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;

@CommandLineProgramProperties(
        summary = "Sanity checks GRIDSS evidence scores by ensuring that both sides of each breakpoint-supporting evidence have the same score",  
        oneLineSummary = "Sanity checks GRIDSS evidence scores.",
        programGroup = gridss.cmdline.programgroups.VariantCalling.class
)
public class SanityCheckEvidence extends FullEvidenceCommandLineProgram {
	private static final Log log = Log.getInstance(SanityCheckEvidence.class);
	@Argument(doc="Margin of error allowed for matching.", optional=true)
	public double ERROR_MARGIN = 0.0001;

	public SanityCheckEvidence() {
		super(false);
	}

	public static void main(String[] argv) {
        System.exit(new SanityCheckEvidence().instanceMain(argv));
    }

	@Override
	protected GridssConfiguration getGridssConfiguration() {
		GridssConfiguration config = super.getGridssConfiguration();
		config.hashEvidenceID = false;
		return config;
	}

	public int sanityCheck(SAMEvidenceSource source) {
		if (source == null) return 0;
		if (!source.getFile().exists()) {
			log.info("Ignoring " + source.getFile());
		}
		try (PairedEvidenceTracker<DirectedEvidence> it = new PairedEvidenceTracker<>("sanitycheck", source.iterator(SAMEvidenceSource.EvidenceSortOrder.EvidenceStartPosition))) {
			while (it.hasNext()) {
				it.next();
			}
			return it.errorCount();
		}
	}
	@Override
	public int doWork(ExecutorService threadpool) {
		for (SAMEvidenceSource ses : getSamEvidenceSources()) {
			log.info("Sanity checking " + ses.getFile().getName());
			if (sanityCheck(ses) != 0) {
				return 1;
			}
		}
		if (sanityCheck(getAssemblySource()) != 0) {
			return 1;
		}
		return 0;
	}
}
