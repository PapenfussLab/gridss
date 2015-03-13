package au.edu.wehi.idsv.sim;

import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.NavigableSet;
import java.util.Random;
import java.util.TreeSet;

import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.IdsvVariantContext;
import au.edu.wehi.idsv.IdsvVariantContextBuilder;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.collect.Lists;
import com.google.common.io.Files;

/**
 * 
 * Simple model of chromothripsis
 * 
 * @author cameron.d
 *
 */
public class FragmentedChromosome extends SimulatedChromosome {
	private final Random rng;
	/**
	 * Linear 1-based genomic coordinate of base immediately before break 
	 */
	private final NavigableSet<Integer> breaks = new TreeSet<Integer>();
	/**
	 * @param reference reference genome
	 * @param breakCleanMargin number of unambiguous bases around the breakpoint
	 */
	public FragmentedChromosome(ProcessingContext context, String chr, int breakCleanMargin, int seed) {
		super(context, chr, breakCleanMargin);
		this.breaks.add(0);
		this.breaks.add(seq.length - 1);
		this.rng = new Random(seed);
	}
	public void breakAt(int position) {
		this.breaks.add(position);
	}
	/**
	 * Shatter the chromosome
	 */
	public void shatter(int pieces) {
		for (int i = 0; i < pieces - 1; i++) {
			randomCleanBreak();
		}
	}
	private int randomCleanBreak() {
		int pos;
		do {
			pos = rng.nextInt(seq.length);
		} while (!isCleanBreak(pos));
		breakAt(pos);
		return pos;
	}
	public void assemble(File fasta, File vcf, int fragments, boolean includeReference) throws IOException {
		List<Fragment> fragList = new ArrayList<Fragment>();
		List<Integer> breakList = Lists.newArrayList(breaks);
		for (int i = 0; i < breaks.size() - 1; i++) {
			int low = breakList.get(i) + 1;
			int high = breakList.get(i + 1);
			//  NNNN
			// 01234
			// ^ ^
			//  ^^
			byte[] fragmentSeq = Arrays.copyOfRange(seq, low - 1, high + 1);
			fragList.add(new Fragment(referenceIndex, low, fragmentSeq, rng.nextBoolean()));
		}
		if (fragments > fragList.size()) throw new IllegalArgumentException("Too many fragments requested");
		Collections.shuffle(fragList, rng);
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
		for (int i = 0; i < fragments; i++) {
			Fragment f = fragList.get(i);
			sb.append(f.getSequence());
			if (last != null) {
				BreakpointSummary bp = new BreakpointSummary(last.getEndBreakend(), f.getStartBreakend());
				String event = String.format("truth_%d_", i);
				calls.add(create(bp, event));
				calls.add(create(bp.remoteBreakpoint(), event));
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
	public IdsvVariantContext create(BreakpointSummary bp, String event) {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(context);
		builder.breakpoint(bp, "")
			.id(event + (bp.isLowBreakend() ? "o" : "h"))
			.attribute(VcfSvConstants.BREAKEND_EVENT_ID_KEY, event)
			.attribute(VcfSvConstants.MATE_BREAKEND_ID_KEY, event + (bp.isLowBreakend() ? "h" : "o"));
		return builder.make();
		
	}
	/**
	 * Determines whether the given breakpoint position contains the required margin of unambiguous, unbroken bases
	 * @param position
	 * @return
	 */
	public boolean isCleanBreak(int position) {
		Integer before = breaks.floor(position);
		Integer after = breaks.ceiling(position);
		if (position - before <= margin || after - position <= margin) return false;
		for (int i = position - margin + 1; i < position + margin; i++) {
			if (!SequenceUtil.isValidBase(seq[i])) {
				return false;
			}
		}
		return true;
	}
}
