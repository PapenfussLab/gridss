package gridss;

import au.edu.wehi.idsv.IdsvVariantContext;
import au.edu.wehi.idsv.VcfBreakendSummary;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import joptsimple.internal.Strings;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;

@CommandLineProgramProperties(
		summary = "Exports single breakend and breakpoint inserted sequences to fasta.",
        oneLineSummary = "Exports single breakend and breakpoint inserted sequences to fasta.",
        programGroup = gridss.cmdline.programgroups.DataConversion.class
)
public class InsertedSequencesToFasta extends CommandLineProgram {
	private static final Log log = Log.getInstance(InsertedSequencesToFasta.class);
	@Argument(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input VCF", optional=false)
	public File INPUT;
	@Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output Fasta", optional=false)
	public File OUTPUT;
	@Argument(doc="Minimum length of inserted sequence to.", optional = true)
	public int MIN_SEQUENCE_LENGTH = 1;

	@Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		try (VCFFileReader vcfReader = new VCFFileReader(INPUT, false)) {
			SAMSequenceDictionary dict = vcfReader.getFileHeader().getSequenceDictionary();
			if (dict == null) {
				throw new RuntimeException("VCF missing sequence headers.");
			}
			try (CloseableIterator<VariantContext> it = vcfReader.iterator()) {
				try (FileOutputStream writer = new FileOutputStream(OUTPUT)) {
					while (it.hasNext()) {
						IdsvVariantContext vc = IdsvVariantContext.create(dict, null, it.next());
						VcfBreakendSummary vbs = new VcfBreakendSummary(dict, vc);
						if (vbs != null && !Strings.isNullOrEmpty(vbs.breakpointSequence) && vbs.breakpointSequence.length() >= MIN_SEQUENCE_LENGTH) {
							StringBuilder sb = new StringBuilder();
							sb.append('>');
							sb.append(vc.getID());
							sb.append('\n');
							sb.append(vbs.breakpointSequence);
							sb.append('\n');
							writer.write(sb.toString().getBytes(StandardCharsets.UTF_8));
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

	public static void main(String[] argv) {
		System.exit(new InsertedSequencesToFasta().instanceMain(argv));
	}
}
