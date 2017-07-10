package gridss;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;

import com.google.common.collect.Iterators;

import au.edu.wehi.idsv.alignment.AlignerFactory;
import au.edu.wehi.idsv.alignment.BreakpointHomology;
import au.edu.wehi.idsv.bed.BedpeIterator;
import au.edu.wehi.idsv.bed.BedpeRecord;
import au.edu.wehi.idsv.bed.BedpeWriter;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

public class AnnotateInexactHomologyBedpe extends ReferenceCommandLineProgram {
	private static final Log log = Log.getInstance(AnnotateInexactHomologyBedpe.class);
	@Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input BEDPE")
    public File INPUT;
	@Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="BEDPE with the score column populated by the length of inexact homology between breakends")
    public File OUTPUT;
	@Option(shortName="D", doc="Number of bases from nominal breakpoint position to consider when calculating homology", optional=true)
	public int DISTANCE = 300;
	@Option(shortName="M", doc="Additional reference bases to include in alignment (to account for indels near the breakpoint).", optional=true)
	public int MARGIN = 32;
	@Option(shortName="UC", doc="1-based index of column containing untemplated sequenced based included in the breakpoint.", optional=true)
	public Integer UNTEMPLATED_SEQUENCE_COLUMN = null;
	protected int doWork() {
		ensureArgs();
		try {
			log.info("Loading aligner");
			AlignerFactory.create();
			annotate();
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
	private void annotate() throws FileNotFoundException, IOException {
		try (BedpeIterator bit = new BedpeIterator(INPUT, getReference())) {
			try (BedpeWriter writer = new BedpeWriter(getReference().getSequenceDictionary(), OUTPUT)) {
				Iterator<InexactHomologyBedpeRecord> it = Iterators.transform(bit, rec -> new InexactHomologyBedpeRecord(rec));
				try (AsyncBufferedIterator<InexactHomologyBedpeRecord> asyncit = new AsyncBufferedIterator<InexactHomologyBedpeRecord>(it, "AsyncAnnotateInexactHomologyBedpe")) {
					while (asyncit.hasNext()) {
						InexactHomologyBedpeRecord rec = asyncit.next();
						writer.write(
								rec.record.bp,
								rec.record.name,
								rec.bh.getLocalHomologyLength() + rec.bh.getRemoteHomologyLength()
								);//Arrays.copyOfRange(rec.record, 10, rec.record.length));
					}
				}
			}
		}
	}
	public static void main(String[] argv) {
        System.exit(new AnnotateInexactHomologyBedpe().instanceMain(argv));
    }
}
