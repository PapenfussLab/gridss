package gridss;

import htsjdk.samtools.*;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

@CommandLineProgramProperties(
		summary = "Exports reads and read pairs with the given names to fastq",
        oneLineSummary = "Exports reads and read pairs with the given names to fastq",
        programGroup = gridss.cmdline.programgroups.DataConversion.class
)
public class ExtractFragmentsToFastq extends CommandLineProgram {
	private static final Log log = Log.getInstance(ExtractFragmentsToFastq.class);
	@Argument(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input file to extract reads.")
	public File INPUT;
	@Argument(doc="File containing read names of matching reads. One per line.")
	public File READ_NAMES;
	@Argument(doc="File to extract unpaired reads to.")
	public File OUTPUT_FQ;
	@Argument(doc="File to extract first read in pair to.")
	public File OUTPUT_FQ1;
	@Argument(doc="File to extract second read in pair to.")
	public File OUTPUT_FQ2;

	public static void main(String[] argv) {
		System.exit(new ExtractFragmentsToFastq().instanceMain(argv));
	}

	private static FastqRecord samToFastq(SAMRecord r) {
		byte[] bases = r.getReadBases().clone();
		byte[] quals = r.getBaseQualities().clone();
		if (r.getReadNegativeStrandFlag()) {
			SequenceUtil.reverseComplement(bases);
			SequenceUtil.reverseQualities(quals);
		}
		return new FastqRecord(
				r.getReadName(),
				new String(bases, StandardCharsets.UTF_8),
				null,
				SAMUtils.phredToFastq(quals));
	}

	@Override
	protected int doWork() {
		Map<String, SAMRecord> lookup = new HashMap<>();
		try {
			HashSet<String> readNames = new HashSet<>(Files.readAllLines(READ_NAMES.toPath()));
			FastqWriterFactory factory = new FastqWriterFactory();
			try (FastqWriter fq1 = factory.newWriter(OUTPUT_FQ1)) {
				try (FastqWriter fq2 = factory.newWriter(OUTPUT_FQ2)) {
					try (FastqWriter fq = factory.newWriter(OUTPUT_FQ)) {
						try (SamReader reader = SamReaderFactory.makeDefault().open(INPUT)) {
							try (SAMRecordIterator it = reader.iterator()) {
								while (it.hasNext()) {
									if (readNames.isEmpty()) {
										log.debug("Found all reads. Stopping input file traversal.");
										break;
									}
									SAMRecord r = it.next();
									String name = r.getReadName();
									if (!r.getSupplementaryAlignmentFlag() && !r.isSecondaryAlignment() && readNames.contains(name)) {
										if (!r.getReadPairedFlag()) {
											fq.write(samToFastq(r));
											readNames.remove(name);
										} else {
											SAMRecord lookupMatch = lookup.get(name);
											if (lookupMatch == null) {
												lookup.put(name, r);
												continue;
											}
											if (lookupMatch.getFirstOfPairFlag() == r.getFirstOfPairFlag()) {
												log.error("Found multiple primary alignment records for %s", name, ". Ignoring all but first.");
												continue;
											}
											SAMRecord r1 = r.getFirstOfPairFlag() ? r : lookupMatch;
											SAMRecord r2 = r.getFirstOfPairFlag() ? lookupMatch : r;
											lookup.remove(name);
											readNames.remove(name);
											fq1.write(samToFastq(r1));
											fq2.write(samToFastq(r2));
										}
									}
								}
							}
						}
						if (!lookup.isEmpty()) {
							log.error("Missing paired primary alignment for ", lookup.size(), " reads. Writing to unpaired fastq.");
							for (SAMRecord r : lookup.values()) {
								fq.write(samToFastq(r));
								readNames.remove(r.getReadName());
							}
						}
						if (!readNames.isEmpty()) {
							log.warn("Missing SAM records for ", readNames.size(), " reads.");
						}
					}
				}
			}
		} catch (IOException e) {
			log.error(e);
			return -1;
		}
		return 0;
	}
}
