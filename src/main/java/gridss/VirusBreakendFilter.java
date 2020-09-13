package gridss;

import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import joptsimple.internal.Strings;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Filters an annotated VIRUSBreakend VCF to only likely integration sites",
        oneLineSummary = "Filters an annotated VIRUSBreakend VCF to only likely integration sites",
        programGroup = gridss.cmdline.programgroups.VariantCalling.class
)
public class VirusBreakendFilter extends CommandLineProgram {
	private static final Log log = Log.getInstance(VirusBreakendFilter.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "VIRUSBreakend VCF file to filter")
	public File INPUT;
	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Filtered VCF")
	public File OUTPUT;
	@Argument(doc = "Minimum portion of host alignment that does not match a simple or low complexity repeat")
	public double MINIMUM_REPEAT_OVERLAP = 0.5;
	@Argument(doc = "Minimum portion of breakend sequence that aligns to host genome")
	public double MINIMUM_HOST_OVERLAP = 0.5;

	public static void main(String[] argv) {
        System.exit(new VirusBreakendFilter().instanceMain(argv));
    }

	@Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		try (VCFFileReader vcfReader = new VCFFileReader(INPUT, false)) {
			VCFHeader header = vcfReader.getFileHeader();
			try (CloseableIterator<VariantContext> it = vcfReader.iterator()) {
				VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
						.setReferenceDictionary(header.getSequenceDictionary())
						.setOutputFile(OUTPUT);
				try (VariantContextWriter vcfWriter = builder.build()) {
					vcfWriter.writeHeader(header);
					while (it.hasNext()) {
						VariantContext vc = it.next();
						if (shouldKeep(vc)) {
							vcfWriter.add(vc);
						}
					}
				}
			}
		}
		return 0;
	}

	private boolean shouldKeep(VariantContext vc) {
		String alt = vc.getAlternateAllele(0).getDisplayString();
		//if (alt.contains("kraken")) return true;
		// single breakends only
		if (!alt.contains(".")) return false;
		int breakendLength = alt.length() - 2;
		List<ChimericAlignment> host = infoToChimeric(vc, VcfInfoAttributes.BREAKEND_ALIGNMENTS.attribute());
		RangeSet<Integer> repeats = repeatRanges(vc);
		for (ChimericAlignment hca : host) {
			int hostAlignmentLength = hca.getLastAlignedBaseReadOffset() - hca.getFirstAlignedBaseReadOffset() + 1;
			if ((hostAlignmentLength / (double)breakendLength) < MINIMUM_HOST_OVERLAP) continue;
			RangeSet<Integer> overlappingRepeats = repeats.subRangeSet(Range.openClosed(hca.getFirstAlignedBaseReadOffset(), hca.getLastAlignedBaseReadOffset() + 1));
			int overlap = overlappingRepeats.asRanges().stream().mapToInt(r -> r.upperEndpoint() - r.lowerEndpoint()).sum();
			if (overlap / (double)hostAlignmentLength < MINIMUM_REPEAT_OVERLAP) return true;
		}
		return false;
	}
	private static List<ChimericAlignment> infoToChimeric(VariantContext vc, String tag) {
		if (Strings.isNullOrEmpty(tag)) return Collections.emptyList();
		List<String> list = vc.getAttributeAsStringList(tag, "");
		if (list.size() == 0 || Strings.isNullOrEmpty(list.get(0))) Collections.emptyList();
		List<ChimericAlignment> calist = list.stream()
				.map(aln -> new ChimericAlignment(aln, "[|:]"))
				.collect(Collectors.toList());
		return calist;
	}
	// half-open
	private static RangeSet<Integer> repeatRanges(VariantContext vc) {
		RangeSet<Integer> rs = TreeRangeSet.create();
		for (ChimericAlignment ca : infoToChimeric(vc, VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_SA_TAG.attribute())) {
			if (ca.rname.contains("Simple_repeat") || ca.rname.contains("Low_complexity")) {
				rs.add(Range.closedOpen(ca.getFirstAlignedBaseReadOffset(), ca.getLastAlignedBaseReadOffset() + 1));
			}
		}
		return rs;
	}
}
