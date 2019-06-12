package gridss;

import au.edu.wehi.idsv.alignment.AlignerFactory;
import au.edu.wehi.idsv.alignment.BreakpointHomology;
import au.edu.wehi.idsv.bed.BedpeIterator;
import au.edu.wehi.idsv.bed.BedpeRecord;
import au.edu.wehi.idsv.bed.BedpeWriter;
import au.edu.wehi.idsv.util.ParallelTransformIterator;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class AnnotateInexactHomologyBedpe extends ReferenceCommandLineProgram {
	private static final Log log = Log.getInstance(AnnotateInexactHomologyBedpe.class);
	@Argument(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input BEDPE")
    public File INPUT;
	@Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="BEDPE with the score column populated by the length of inexact homology between breakends")
    public File OUTPUT;
	@Argument(shortName="D", doc="Number of bases from nominal breakpoint position to consider when calculating homology", optional=true)
	public int DISTANCE = 300;
	@Argument(shortName="M", doc="Additional reference bases to include in alignment (to account for indels near the breakpoint).", optional=true)
	public int MARGIN = 32;
	@Argument(shortName="UC", doc="1-based index of column containing untemplated sequenced based included in the breakpoint.", optional=true)
	public Integer UNTEMPLATED_SEQUENCE_COLUMN = null;
	@Argument(doc="Number of worker threads to spawn. Defaults to number of cores available."
			+ " Note that I/O threads are not included in this worker thread count so CPU usage can be higher than the number of worker thread.",
    		shortName="THREADS")
    public int WORKER_THREADS = Runtime.getRuntime().availableProcessors();
	protected int doWork() {
		ensureArgs();
		try {
			log.info("Loading aligner");
			AlignerFactory.create();
			log.info(String.format("Using %d worker threads", WORKER_THREADS));
    		ExecutorService threadpool = Executors.newFixedThreadPool(WORKER_THREADS, new ThreadFactoryBuilder().setDaemon(false).setNameFormat("Worker-%d").build());
			annotate(threadpool);
			threadpool.shutdown();
		} catch (Exception e) {
			log.error(e);
			return 1;
		}
		return 0;
	}
	private void ensureArgs() {
		IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
	}
	public class InexactHomologyBedpeRecord {
		public final BedpeRecord record;
		public final BreakpointHomology bh;
		public InexactHomologyBedpeRecord(BedpeRecord record) {
			this.record = record;
			String untemplated = "";
			if (UNTEMPLATED_SEQUENCE_COLUMN != null && record.fields.length >= UNTEMPLATED_SEQUENCE_COLUMN) {
				untemplated = record.fields[UNTEMPLATED_SEQUENCE_COLUMN - 1];
			}
			bh = BreakpointHomology.calculate(getReference(), record.bp.getNominalPosition(), untemplated, DISTANCE, MARGIN);
		}
	}
	private void annotate(ExecutorService threadpool) throws FileNotFoundException, IOException {
		SAMSequenceDictionary dict = getReference().getSequenceDictionary();
		try (BedpeIterator bit = new BedpeIterator(INPUT, dict)) {
			try (BedpeWriter writer = new BedpeWriter(dict, OUTPUT)) {
				ParallelTransformIterator<BedpeRecord, InexactHomologyBedpeRecord> asyncit = new ParallelTransformIterator<BedpeRecord, InexactHomologyBedpeRecord>(
						bit, rec -> new InexactHomologyBedpeRecord(rec), WORKER_THREADS + 1, threadpool);
				while (asyncit.hasNext()) {
					InexactHomologyBedpeRecord rec = asyncit.next();
					writer.write(
							rec.record.bp,
							rec.record.name,
							Integer.toString(rec.bh.getLocalHomologyLength() + rec.bh.getRemoteHomologyLength()),
							new String[] {
									Integer.toString(rec.bh.getLocalHomologyLength()),
									Integer.toString(rec.bh.getRemoteHomologyLength()),
							}
							);//Arrays.copyOfRange(rec.record, 10, rec.record.length));
				}
			}
		}
	}
	public static void main(String[] argv) {
        System.exit(new AnnotateInexactHomologyBedpe().instanceMain(argv));
    }
}
