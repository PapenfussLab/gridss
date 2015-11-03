package au.edu.wehi.idsv.sim;

import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.sim.SequentialVariantPlacer.ContigExhaustedException;
import au.edu.wehi.idsv.vcf.SvType;
import au.edu.wehi.idsv.vcf.VcfFilter;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.io.Files;

public class SimpleVariantChromosome extends SimulatedChromosome {
	private static class Event {
		public Event(SvType t, int s, int z) {
			type = t;
			start = s;
			size = z;
		}
		public final SvType type;
		/**
		 * Genomic position immediately prior to event
		 */
		public final int start;
		/**
		 * Event size in base pairs
		 */
		public final int size;
	}
	private RandomBaseGenerator baseGen;
	/**
	 * @param reference reference genome
	 * @param margin number of unambiguous bases around the breakpoint
	 */
	public SimpleVariantChromosome(ProcessingContext context, String chr, int margin, int seed) {
		super(context, chr, margin, seed);
		this.baseGen = new RandomBaseGenerator(seed);
	}
	public void assemble(File fasta, File vcf, boolean includeReference, List<SvType> type, List<Integer> size, int countPerEventTypeSize) throws IOException {
		List<Event> variantList = populateEvents(type, size, countPerEventTypeSize);
		// ensure counts are the same
		variantList.subList(0, variantList.size() - variantList.size() % countPerEventTypeSize);
		List<VariantContext> list = new ArrayList<VariantContext>();
		
		StringBuilder sb = new StringBuilder(">variant." + getChr() + "\n");
		
		int genomicPosition = 0; // emitting up to and including this genomic position
		for (Event e : variantList) {
			//  012345678   array index
			// 0123456789   genomic position
			// ^ ||  |		genomicPosition = last emitted position 
			//   ^|  |		start
			//    ^^^^		event bases
			// emit up to start location exclusive
			String beforeAnchorSequence = new String(seq, genomicPosition, e.start - genomicPosition - 1);
			sb.append(beforeAnchorSequence);
			genomicPosition = e.start;
			String refBefore =  new String(seq, genomicPosition - 1, 1);
			String altSeq = getVariantSeq(e.type, e.size, genomicPosition, 1, 0);
			sb.append(altSeq);
			VariantContextBuilder builder = new VariantContextBuilder();
			builder.id(String.format("%s.%d.%s%d", getChr(), e.start, e.type, e.size))
				.chr(getChr())
				.start(e.start)
				.stop(e.start)
				.attribute(VcfSvConstants.SV_LENGTH_KEY, e.type == SvType.DEL ? -e.size : e.size)
				.attribute(VcfSvConstants.END_KEY, genomicPosition)
				.attribute(VcfSvConstants.SV_TYPE_KEY, e.type)
				//.attribute("ALT", altSeq)
				.alleles(refBefore, "<" + e.type + ">");
			if (altSeq.equalsIgnoreCase(new String(seq, genomicPosition - 1, getGenomicWidth(e.type, e.size) + 1))) {
				builder.filter(VcfFilter.REFERENCE_ALLELE.filter());
			} else {
				int homBefore = homLenBefore(e.type, e.size, e.start);
				int homAfter = homLenAfter(e.type, e.size, e.start);
				if (homBefore > 0 || homAfter > 0) {
					builder.attribute(VcfSvConstants.HOMOLOGY_LENGTH_KEY, homBefore + homAfter);
					builder.attribute(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY, new int[] { -homBefore, homAfter });
				}
			}
			list.add(builder.make());
			genomicPosition += getGenomicWidth(e.type, e.size);
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
		Files.write(sb.toString(), fasta, StandardCharsets.US_ASCII);
		writeVcf(vcf, list);
	}
	private int homLenBefore(SvType t, int size, int genomicPosition) {
		int len = 0;
		while (genomicPosition - (len + 1) > 0 &&
				getVariantSeq(t, size, genomicPosition, len + 1, 0).equals(
				getVariantSeq(t, size, genomicPosition - (len + 1), 0, len + 1))) {
			len++;
		}
		return len;
	}
	private int homLenAfter(SvType t, int size, int genomicPosition) {
		int len = 0;
		while (genomicPosition + (len + 1) < seq.length - 1 &&
				getVariantSeq(t, size, genomicPosition, 0, len + 1).equals(
				getVariantSeq(t, size, genomicPosition + (len + 1), len + 1, 0))) {
			len++;
		}
		return len;
	}
	private String getVariantSeq(SvType t, int size, int genomicPosition, int basesBefore, int basesAfter) {
		String ref = new String(seq, genomicPosition, size);
		String altSeq;
		switch (t) {
		case INS:
			altSeq = new String(baseGen.getBases(size));
			break;
		case DEL:
			altSeq = "";
			break;
		case INV:
			altSeq = SequenceUtil.reverseComplement(ref);
			break;
		case DUP:
			altSeq = ref + ref;
			break;
		default:
			throw new RuntimeException("NYI");
		}
		String before = new String(seq, genomicPosition - basesBefore, basesBefore);
		String after = new String(seq, genomicPosition + getGenomicWidth(t, size), basesAfter);
		return before + altSeq + after;
	}
	private List<Event> populateEvents(List<SvType> typeList, List<Integer> sizeList, int max) {
		SequentialVariantPlacer placer = new SequentialVariantPlacer(seq, margin);
		List<Event> variantList = new ArrayList<Event>();
		try {
			for (int i = 0; i < max; i++) {
				for (int size : sizeList) {
					for (SvType type : typeList) {
						int width = getGenomicWidth(type, size);
						variantList.add(new Event(type, placer.getNext(width + 1), size));
					}
				}
			}
		} catch (ContigExhaustedException e) {
		}
		return variantList;
	}
	private int getGenomicWidth(SvType type, int size) {
		switch (type) {
		case INS:
			return 0;
		case DEL:
		case INV:
		case DUP:
			return size;
		case BND:
		case CNV:
		default:
			throw new RuntimeException("Not implemented by this simulator");
		}
	}
}
