package au.edu.wehi.idsv.sim;

import au.edu.wehi.idsv.GenomicProcessingContext;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 * 
 * Simple model of chromothripsis
 * 
 * @author Daniel Cameron
 *
 */
public class FragmentedChromosome extends SimulatedChromosome {
	private static final Log log = Log.getInstance(FragmentedChromosome.class);
	protected final int minFragmentLength;
	protected final int maxFragmentLength;
	private final int telomereLength;

	/**
	 * @param breakMargin number of unambiguous bases around the breakpoint
	 */
	public FragmentedChromosome(GenomicProcessingContext context, String chr, int telomereLength, int breakMargin, int minFragmentLength, int maxFragmentLength, int seed) {
		super(context, chr, breakMargin, seed);
		this.minFragmentLength = minFragmentLength;
		this.maxFragmentLength = maxFragmentLength;
		this.telomereLength = telomereLength;
	}
	private class RandomFragmentIterator implements Iterator<Fragment> {
		@Override
		public boolean hasNext() {
			return true;
		}
		@Override
		public Fragment next() {
			int fragmentLength = getRandomFragmentLength();
			return createFragment(1 + rng.nextInt(seq.length - fragmentLength), fragmentLength, rng.nextBoolean());
		}
		@Override
		public void remove() {
		}
	}
	public int getRandomFragmentLength() {
		int fragmentLength = minFragmentLength + rng.nextInt(maxFragmentLength - minFragmentLength + 1);
		return fragmentLength;
	}
	protected Iterator<Fragment> candidateFragments() {
		return new RandomFragmentIterator();
	}
	public void assemble(File fasta, File vcf, int fragments, boolean includeReference) throws IOException {
		List<Fragment> fragList = new ArrayList<Fragment>();
		RangeSet<Integer> invalid = calcInvalidBreakPositions(margin + maxFragmentLength);
		Iterator<Fragment> it = candidateFragments();
		while (it.hasNext() && fragList.size() < fragments) {
			Fragment f = it.next();
			if (f == null || !invalid.subRangeSet(Range.closedOpen(f.getLowBreakend().start, f.getHighBreakend().end)).isEmpty()) {
				// skip fragments we can't use
				continue;
			}
			invalid.add(Range.closedOpen(f.getLowBreakend().start - margin, f.getHighBreakend().end + 1 + margin));
			fragList.add(f);
		}
		if (telomereLength > 0) {
			fragList.add(0, createFragment(1, telomereLength, false));
			fragList.add(createFragment(seq.length - telomereLength + 1, telomereLength, false));
		}
		log.info(String.format("%d fragments created", fragList.size()));
		assemble(fasta, vcf, fragList, includeReference);
	}
}
