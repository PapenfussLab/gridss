package au.edu.wehi.idsv.sim;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.bed.BedpeRecord;
import au.edu.wehi.idsv.bed.BedpeWriter;
import au.edu.wehi.idsv.picard.BufferedReferenceSequenceFile;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.vcf.*;
import com.google.common.collect.Lists;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;
import com.google.common.io.Files;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.*;

public class SimulatedChromosome {
	protected final Random rng;
	protected final byte[] seq;
	protected final int referenceIndex;
	protected final int margin;
	protected ReferenceLookup ref;
	protected final String chr;
	protected final GenomicProcessingContext context;
	protected Fragment createFragment(int genomicStart, int length, boolean reversed) {
		return new Fragment(referenceIndex, genomicStart, Arrays.copyOfRange(seq, genomicStart - 1,  genomicStart - 1 + length), reversed);
	}
	/**
	 * @param context reference genome
	 * @param margin number of unambiguous bases around the breakpoint
	 */
	public SimulatedChromosome(GenomicProcessingContext context, String chr, int margin, int seed) {
		this.context = context;
		this.ref = new BufferedReferenceSequenceFile(context.getReference());
		this.referenceIndex = ref.getSequenceDictionary().getSequence(chr).getSequenceIndex();
		this.chr = ref.getSequenceDictionary().getSequence(chr).getSequenceName();
		this.seq = ref.getSequence(chr).getBases();
		this.margin = margin;
		this.rng = new Random(seed);
	}
	protected void writeVcf(File vcf, Iterable<VariantContext> calls) {
		VariantContextWriterBuilder builder = context.getVariantContextWriterBuilder(vcf, true);
		VariantContextWriter writer = builder.build();
		VCFHeader header = new VCFHeader();
		// Standard SV headers we use
		SimpleEvent.addVcfHeaders(header);
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.CONFIDENCE_INTERVAL_START_POSITION);
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.HOMOLOGY_LENGTH);
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.HOMOLOGY_SEQUENCE);
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.MATE_BREAKEND_ID);
		header.addMetaDataLine(VcfStructuralVariantHeaderLines.EVENT_ID);
		writer.writeHeader(header);
		for (VariantContext vc : calls) {
			writer.add(vc);
		}
		writer.close();
	}
	protected String getChr() { return chr; }
	/**
	 * Calculates the ranges where a breakpoint would not have sufficient margin on either side
	 * @param basesBuffer
	 * @return invalid breakpoints
	 */
	protected RangeSet<Integer> calcInvalidBreakPositions(int basesBuffer) {
		RangeSet<Integer> invalid = TreeRangeSet.create();
		int firstAmbiguous = 1;
		boolean inRange = false;
		for (int i = 0; i < seq.length; i++) {
			int genomicPos = i + 1;
			if (!SequenceUtil.isValidBase(seq[i])) {
				if (!inRange) {
					firstAmbiguous = genomicPos;
					inRange = true;
				}
			} else if (inRange) {
				inRange = false;
				invalid.add(Range.closedOpen(firstAmbiguous - basesBuffer, genomicPos + basesBuffer));
			}
		}
		if (inRange) {
			invalid.add(Range.closedOpen(firstAmbiguous - basesBuffer, seq.length + basesBuffer));
		}
		invalid.add(Range.closedOpen(0, basesBuffer + 1));
		invalid.add(Range.closedOpen(seq.length - basesBuffer - 1, seq.length));
		return invalid;
	}
	protected SAMSequenceDictionary assemble(File fasta, File vcf, File bedpe, List<Fragment> fragList, boolean includeReference) throws IOException {
		List<IdsvVariantContext> calls = Lists.newArrayList();
		List<BedpeRecord> breakpointLocations = new ArrayList<>();
		Fragment last = null;
		int lastPosition = 0;
		StringBuilder variantSeq = new StringBuilder();
		for (int i = 0; i < fragList.size(); i++) {
			Fragment f = fragList.get(i);
			variantSeq.append(f.getSequence());
			int position = variantSeq.length();
			if (last != null) {
				BreakpointSummary bp = new BreakpointSummary(last.getEndBreakend(), f.getStartBreakend());
				String event = String.format("truth_%d_", i);
				calls.add(create(bp, event).make());
				calls.add(create(bp.remoteBreakpoint(), event).make());
				breakpointLocations.add(new BedpeRecord(event,
						new BreakpointSummary(0, BreakendDirection.Forward, lastPosition, 0, BreakendDirection.Backward, lastPosition + 1),
						"1"));
			}
			last = f;
			lastPosition = position;
		}
		SAMSequenceDictionary dict = new SAMSequenceDictionary();
		String variantContigName = "chromothripsis." + getChr();
		StringBuilder sb = new StringBuilder();
		sb.append(">" + variantContigName + "\n");
		sb.append(variantSeq);
		dict.addSequence(new SAMSequenceRecord(variantContigName, variantSeq.length()));
		if (includeReference) {
			sb.append(">");
			sb.append(getChr());
			sb.append("\n");
			sb.append(new String(seq, StandardCharsets.US_ASCII));
			sb.append("\n");
			dict.addSequence(new SAMSequenceRecord(getChr(), seq.length));
		}
		Collections.sort(calls, IdsvVariantContext.ByLocationStart);
		Files.asCharSink(fasta, StandardCharsets.US_ASCII).write(sb.toString());
		try (VariantContextWriter writer = context.getVariantContextWriter(vcf, true)) {
			for (VariantContext vc : calls) {
				writer.add(vc);
			}
		}
		if (bedpe != null) {
			try (BedpeWriter writer = new BedpeWriter(dict, bedpe)) {
				writer.writeHeader(false, false);
				for (BedpeRecord r : breakpointLocations) {
					writer.write(r);
				}
			}
		}
		return dict;
	}
	protected IdsvVariantContextBuilder create(BreakpointSummary bp, String event) {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(context);
		builder.breakpoint(bp, "")
			.id(event + (bp.isLowBreakend() ? "o" : "h"))
			.attribute(VcfSvConstants.EVENT_ID_KEY, event)
			.attribute(VcfSvConstants.MATE_BREAKEND_ID_KEY, event + (bp.isLowBreakend() ? "h" : "o"))
			.attribute(VcfSvConstants.MATE_BREAKEND_ID_KEY, event + (bp.isLowBreakend() ? "h" : "o"));
		return builder;
	}
}