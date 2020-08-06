package gridss.kraken;

import au.edu.wehi.idsv.kraken.AnnotateKraken;
import au.edu.wehi.idsv.kraken.KrakenParser;
import au.edu.wehi.idsv.vcf.GridssVcfConstants;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.*;

@CommandLineProgramProperties(
        summary = "Annotates breakpoint variant calls with Kraken2 classifications",
        oneLineSummary = "Annotates breakpoint variant calls Kraken2 classifications",
        programGroup = gridss.cmdline.programgroups.VariantCalling.class
)
public class AnnotateVariantsKraken extends CommandLineProgram {
	private static final Log log = Log.getInstance(AnnotateVariantsKraken.class);
	@Argument(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input VCF file")
	public File INPUT;
	@Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output VCF file")
	public File OUTPUT;
	@Argument(shortName = "K", doc="Kraken2 output file. Records must be in the same order as the VCF.")
	public File KRAKEN_INPUT;

	@Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsReadable(KRAKEN_INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		try (KrakenParser parser = new KrakenParser(new BufferedReader(new InputStreamReader(new FileInputStream(INPUT))))) {
			try (VCFFileReader vcfReader = new VCFFileReader(INPUT, false)) {
				VCFHeader header = vcfReader.getFileHeader();
				if (header.getSequenceDictionary() == null) {
					throw new RuntimeException("INPUT VCF missing sequence definitions.");
				}
				try (CloseableIterator<VariantContext> it = vcfReader.iterator()) {
					VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
							.setReferenceDictionary(header.getSequenceDictionary())
							.setOutputFile(OUTPUT);
					try (VariantContextWriter vcfWriter = builder.build()) {
						header.addMetaDataLine(VcfInfoAttributes.INSERTED_SEQUENCE_NCBI_TAXONOMY_ID.infoHeader());
						vcfWriter.writeHeader(header);
						try (AnnotateKraken ak = new AnnotateKraken(new KrakenParser(new BufferedReader(new InputStreamReader(new FileInputStream(KRAKEN_INPUT)))), it)) {
							while (ak.hasNext()) {
								VariantContext vc = ak.next();
								vcfWriter.add(vc);
							}
						}
					}
				}
			}
		} catch (IOException e) {
			log.error(e);
			throw new RuntimeIOException(e);
		}
		return 0;
	}

	public static void main(String[] argv) {
		System.exit(new AnnotateVariantsKraken().instanceMain(argv));
	}
}
