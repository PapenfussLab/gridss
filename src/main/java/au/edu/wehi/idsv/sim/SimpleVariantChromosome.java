package au.edu.wehi.idsv.sim;

import au.edu.wehi.idsv.GenomicProcessingContext;
import au.edu.wehi.idsv.sim.SequentialVariantPlacer.ContigExhaustedException;
import au.edu.wehi.idsv.vcf.SvType;
import au.edu.wehi.idsv.vcf.VcfFilter;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import com.google.common.io.Files;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import static au.edu.wehi.idsv.vcf.SvType.INS;

public class SimpleVariantChromosome extends SimulatedChromosome {
	private final boolean useSymbolicAllele;
	private final RandomBaseGenerator baseGen;
	/**
	 * @param margin number of unambiguous bases around the breakpoint
	 */
	public SimpleVariantChromosome(GenomicProcessingContext context, String chr, int margin, int seed, boolean useSymbolicAllele) {
		super(context, chr, margin, seed);
		this.baseGen = new RandomBaseGenerator(seed);
		this.useSymbolicAllele = useSymbolicAllele;
	}
	public void assemble(File fasta, File vcf, boolean includeReference, List<SvType> type, List<Integer> size, int countPerEventTypeSize) throws IOException {
		List<SimpleEvent> variantList = populateEvents(type, size, countPerEventTypeSize);
		// ensure counts are the same
		variantList.subList(0, variantList.size() - variantList.size() % countPerEventTypeSize);
		variantList.sort(Comparator.comparingInt(o -> o.start));
		List<VariantContext> list = new ArrayList<VariantContext>();
		StringBuilder sb = new StringBuilder(">variant." + getChr() + "\n");

		int genomicPosition = 0; // emitted up to and including this genomic position
		for (SimpleEvent e : variantList) {
			String beforeVariantSequence = ref.getSubsequenceAt(chr, genomicPosition + 1, e.start).getBaseString();
			String variantSequence = e.getVariantSeq(ref, 0, 0);
			sb.append(beforeVariantSequence);
			sb.append(variantSequence);
			list.add(e.asVariantContextBuilder(ref, useSymbolicAllele).make());
			genomicPosition = e.start + e.getGenomicWidth();
		}
		sb.append(new String(seq, genomicPosition, margin));
		genomicPosition += margin;
		if (includeReference) {
			sb.append("\n>");
			sb.append(getChr());
			sb.append("\n");
			sb.append(new String(seq, 0, genomicPosition));
			sb.append("\n");
		}
		Files.asCharSink(fasta, StandardCharsets.US_ASCII).write(sb.toString());
		writeVcf(vcf, list);
	}
	private List<SimpleEvent> populateEvents(List<SvType> typeList, List<Integer> sizeList, int max) {
		SequentialVariantPlacer placer = new SequentialVariantPlacer(seq, margin);
		List<SimpleEvent> variantList = new ArrayList<SimpleEvent>();
		try {
			for (int i = 0; i < max; i++) {
				for (int size : sizeList) {
					for (SvType type : typeList) {
						int width = SimpleEvent.getGenomicWidth(type, size);
						String insSeq = type != INS ? "" : new String(baseGen.getBases(size));
						variantList.add(new SimpleEvent(type, referenceIndex, size, placer.getNext(width + 1), insSeq));
					}
				}
			}
		} catch (ContigExhaustedException e) {
		}
		return variantList;
	}
}
