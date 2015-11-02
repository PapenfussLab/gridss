package au.edu.wehi.idsv.sim;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;

import au.edu.wehi.idsv.ProcessingContext;
import htsjdk.samtools.util.Log;

/**
 * 
 * Simple model of chromothripsis
 * 
 * @author cameron.d
 *
 */
public class FragmentedChromosome extends SimulatedChromosome {
	private static final Log log = Log.getInstance(FragmentedChromosome.class);
	protected final int fragmentLength;
	/**
	 * @param reference reference genome
	 * @param breakMargin number of unambiguous bases around the breakpoint
	 */
	public FragmentedChromosome(ProcessingContext context, String chr, int breakMargin, int fragmentLength, int seed) {
		super(context, chr, breakMargin, seed);
		this.fragmentLength = fragmentLength;
	}
	private class RandomFragmentIterator implements Iterator<Fragment> {
		@Override
		public boolean hasNext() {
			return true;
		}
		@Override
		public Fragment next() {
			return createFragment(1 + rng.nextInt(seq.length - fragmentLength), fragmentLength, rng.nextBoolean());
		}
		@Override
		public void remove() {
		}
	}
	protected Iterator<Fragment> candidateFragments() {
		return new RandomFragmentIterator();
	}
	public void assemble(File fasta, File vcf, int fragments, boolean includeReference) throws IOException {
		List<Fragment> fragList = new ArrayList<Fragment>();
		RangeSet<Integer> invalid = calcInvalidBreakPositions(margin + fragmentLength);
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
		log.info(String.format("%d fragments created", fragList.size()));
		assemble(fasta, vcf, fragList, includeReference);
	}
}
