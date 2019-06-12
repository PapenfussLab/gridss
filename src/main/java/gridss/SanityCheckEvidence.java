package gridss;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.SAMEvidenceSource;
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
	public void sanityCheck(SAMEvidenceSource source) {
		if (source == null) return;
		if (!source.getFile().exists()) {
			log.info("Ignoring " + source.getFile());
		}
		try (CloseableIterator<DirectedEvidence> it = new PairedEvidenceTracker<DirectedEvidence>("sanitycheck", source.iterator())) {
			while (it.hasNext()) {
				it.next();
			}
		}
	}
	@Override
	public int doWork(ExecutorService threadpool) throws IOException, InterruptedException, ExecutionException {
		for (SAMEvidenceSource ses : getSamEvidenceSources()) {
			sanityCheck(ses);
		}
		sanityCheck(getAssemblySource());
		return 0;
	}
}
