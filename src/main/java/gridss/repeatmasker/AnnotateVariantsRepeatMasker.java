package gridss.repeatmasker;

import au.edu.wehi.idsv.kraken.KrakenParser;
import au.edu.wehi.idsv.repeatmasker.AnnotateRepeatMasker;
import au.edu.wehi.idsv.repeatmasker.RepeatMaskerCodec;
import au.edu.wehi.idsv.repeatmasker.RepeatMaskerFeature;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.readers.LineIterator;
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
import java.util.List;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Annotates single breakend and breakpoint inserted sequences with RepeatMasker classifications." +
				" Since RepeatMasker does not retain record ordering, this program must be given enough memory to" +
				" load the entire RepeatMasker output file into memory.",
        oneLineSummary = "Annotates single breakend and breakpoint inserted sequences with RepeatMasker classifications",
        programGroup = gridss.cmdline.programgroups.VariantCalling.class
)
public class AnnotateVariantsRepeatMasker extends CommandLineProgram {
	private static final Log log = Log.getInstance(AnnotateVariantsRepeatMasker.class);
	@Argument(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input VCF file")
	public File INPUT;
	@Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output VCF file")
	public File OUTPUT;
	@Argument(shortName = "RM", doc="RepeatMasker output or detailed alignment file." +
			" If an alignment file is specified, the INSRM CIGAR and edit distance will be fully populated." +
			" If an output file is specifed, only a minimal CIGAR will be reported as the full alignment is not available")
	public File REPEAT_MASKER;
	@Argument(doc="INFO fields to populate. Valid values are INSRM, INSRMRT, INSRMRC, INSRMRO, INSRMP")
	public List<String> TAGS = AnnotateRepeatMasker.REPEAT_MASKER_ATTRIBUTES.stream().map(x -> x.attribute()).collect(Collectors.toList());
	@Override
	public int doWork() {
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsReadable(REPEAT_MASKER);
		IOUtil.assertFileIsWritable(OUTPUT);
		try (KrakenParser parser = new KrakenParser(new BufferedReader(new InputStreamReader(new FileInputStream(INPUT))))) {
			try (VCFFileReader vcfReader = new VCFFileReader(INPUT, false)) {
				VCFHeader header = vcfReader.getFileHeader();
				if (header.getSequenceDictionary() == null) {
					throw new RuntimeException("INPUT VCF missing sequence definitions.");
				}
				log.info("Loading ", REPEAT_MASKER);
				Multimap<String, RepeatMaskerFeature> lookup = TreeMultimap.create(String::compareTo, RepeatMaskerFeature.ByUniqueID);
				try (AbstractFeatureReader<RepeatMaskerFeature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(REPEAT_MASKER.getPath(), new RepeatMaskerCodec(), false)) {
					try(CloseableTribbleIterator<RepeatMaskerFeature> rmit = reader.iterator()) {
						while (rmit.hasNext()) {
							RepeatMaskerFeature f = rmit.next();
							lookup.put(f.getContig(), f);
					}
				}
				log.info("Loaded ", lookup.size(), " RepeatMasker records.");
				log.info("Annotating ", OUTPUT);
					try (CloseableIterator<VariantContext> it = vcfReader.iterator()) {
						VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
								.setReferenceDictionary(header.getSequenceDictionary())
								.setOutputFile(OUTPUT);
						try (VariantContextWriter vcfWriter = builder.build()) {
							for (VcfInfoAttributes a : AnnotateRepeatMasker.REPEAT_MASKER_ATTRIBUTES) {
								if (TAGS.contains(a.attribute())) {
									header.addMetaDataLine(a.infoHeader());
								}
							}
							vcfWriter.writeHeader(header);
							while (it.hasNext()) {
								VariantContext vc = it.next();
								vc = AnnotateRepeatMasker.annotate(header.getSequenceDictionary(), TAGS, vc, lookup.removeAll(vc.getID()));
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
		System.exit(new AnnotateVariantsRepeatMasker().instanceMain(argv));
	}
}
