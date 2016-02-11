package au.edu.wehi.idsv.sim;

import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.GenomicProcessingContext;
import au.edu.wehi.idsv.IdsvVariantContext;
import au.edu.wehi.idsv.IdsvVariantContextBuilder;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.collect.Lists;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;
import com.google.common.io.Files;

public class SimulatedChromosome {
	protected final Random rng;
	protected final byte[] seq;
	protected final int referenceIndex;
	protected final int margin;
	protected final ReferenceLookup ref;
	protected final String chr;
	protected final GenomicProcessingContext context;
	protected Fragment createFragment(int genomicStart, int length, boolean reversed) {
		return new Fragment(referenceIndex, genomicStart, Arrays.copyOfRange(seq, genomicStart - 1,  genomicStart - 1 + length), reversed);
	}
	/**
	 * @param reference reference genome
	 * @param breakCleanMargin number of unambiguous bases around the breakpoint
	 */
	public SimulatedChromosome(GenomicProcessingContext context, String chr, int margin, int seed) {
		this.context = context;
		this.ref = context.getReference();
		this.referenceIndex = ref.getSequenceDictionary().getSequence(chr).getSequenceIndex();
		this.chr = ref.getSequenceDictionary().getSequence(chr).getSequenceName();
		this.seq = ref.getSequence(chr).getBases();
		this.margin = margin;
		this.rng = new Random(seed);
	}
	protected void writeVcf(File vcf, Iterable<VariantContext> calls) {
		VariantContextWriter writer = context.getVariantContextWriter(vcf, true);
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
	protected void assemble(File fasta, File vcf, List<Fragment> fragList, boolean includeReference) throws IOException {
		StringBuilder sb = new StringBuilder();
		if (includeReference) {
			sb.append(">");
			sb.append( getChr());
			sb.append("\n");
			sb.append(new String(seq, StandardCharsets.US_ASCII));
			sb.append("\n");
		}
		sb.append(">chromothripsis." + getChr() + "\n");
		List<IdsvVariantContext> calls = Lists.newArrayList();
		Fragment last = null;
		for (int i = 0; i < fragList.size(); i++) {
			Fragment f = fragList.get(i);
			sb.append(f.getSequence());
			if (last != null) {
				BreakpointSummary bp = new BreakpointSummary(last.getEndBreakend(), f.getStartBreakend());
				String event = String.format("truth_%d_", i);
				calls.add(create(bp, event).make());
				calls.add(create(bp.remoteBreakpoint(), event).make());
			}
			last = f;
		}
		Collections.sort(calls, IdsvVariantContext.ByLocationStart);
		Files.write(sb.toString(), fasta, StandardCharsets.US_ASCII);
		VariantContextWriter writer = context.getVariantContextWriter(vcf, true);
		for (VariantContext vc : calls) {
			writer.add(vc);
		}
		writer.close();
	}
	protected IdsvVariantContextBuilder create(BreakpointSummary bp, String event) {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(context);
		builder.breakpoint(bp, "")
			.id(event + (bp.isLowBreakend() ? "o" : "h"))
			.attribute(VcfSvConstants.BREAKEND_EVENT_ID_KEY, event)
			.attribute(VcfSvConstants.MATE_BREAKEND_ID_KEY, event + (bp.isLowBreakend() ? "h" : "o"));
		return builder;
	}
}