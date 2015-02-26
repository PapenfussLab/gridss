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

import com.google.common.collect.Lists;
import com.google.common.io.Files;

/**
 * 
 * Simple model of chromothripsis
 * 
 * @author cameron.d
 *
 */
public class FragmentedChromosome {
	private final byte[] seq;
	private final int referenceIndex;
	private final int margin;
	private final Random rng;
	/**
	 * Linear 1-based genomic coordinate of base immediately before break 
	 */
	private final NavigableSet<Integer> breaks = new TreeSet<Integer>();
	private final ProcessingContext context;
	/**
	 * @param reference reference genome
	 * @param breakCleanMargin number of unambiguous bases around the breakpoint
	 */
	public FragmentedChromosome(ProcessingContext context, String chr, int breakCleanMargin, int seed) {
		this.context = context;
		this.referenceIndex = context.getReference().getSequenceDictionary().getSequence(chr).getSequenceIndex();
		this.seq = context.getReference().getSequence(chr).getBases();
		this.margin = breakCleanMargin;
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
	public void assemble(File fasta, File vcf, int fragments) throws IOException {
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
		StringBuilder sb = new StringBuilder(">chromthripsis\n");
		List<IdsvVariantContext> calls = Lists.newArrayList();
		Fragment last = null;
		for (int i = 0; i < fragments; i++) {
			Fragment f = fragList.get(i);
			sb.append(f.getSequence());
			if (last != null) {
				BreakpointSummary bp = new BreakpointSummary(last.getEndBreakend(), f.getStartBreakend());
				calls.add(new IdsvVariantContextBuilder(context).breakpoint(bp, "").make());
				calls.add(new IdsvVariantContextBuilder(context).breakpoint(bp.remoteBreakpoint(), "").make());
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
