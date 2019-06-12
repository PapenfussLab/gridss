package au.edu.wehi.idsv.sim;

import au.edu.wehi.idsv.GenomicProcessingContext;
import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class RepeatFragmentedChromosome extends FragmentedChromosome {
	private static final Log log = Log.getInstance(RepeatFragmentedChromosome.class);
	private File repeatMaskerOutput;
	private String repeatClassFamily;
	public RepeatFragmentedChromosome(GenomicProcessingContext context, String chr, int breakMargin, int fragmentLength, File repeatMaskerOutput, String repeatClassFamily, int seed) {
		super(context, chr, breakMargin, fragmentLength, seed);
		this.repeatMaskerOutput = repeatMaskerOutput;
		this.repeatClassFamily = repeatClassFamily;
	}
	private static RangeMap<Integer, RepeatMaskerRepeat> toLookup(List<RepeatMaskerRepeat> repeatList) throws IOException {
		RangeMap<Integer, RepeatMaskerRepeat> elements = TreeRangeMap.create();
		for (RepeatMaskerRepeat repeat : repeatList) {
			elements.put(Range.closed(repeat.begin, repeat.end), repeat);
		}
		return elements;
	}
	@Override
	protected Iterator<Fragment> candidateFragments() {
		try {
			log.info(String.format("Loading %s from %s", repeatClassFamily, repeatMaskerOutput));
			final List<RepeatMaskerRepeat> repeatList = RepeatMaskerRepeat.loadRepeatMaskerOutput(repeatMaskerOutput, getChr(), repeatClassFamily);
			final RangeMap<Integer, RepeatMaskerRepeat> repeats = toLookup(repeatList);
			Collections.shuffle(repeatList, rng);
			return Iterators.transform(repeatList.iterator(), new Function<RepeatMaskerRepeat, Fragment>() {
				@Override
				public Fragment apply(RepeatMaskerRepeat r) {
					Fragment f = null;
					if (rng.nextBoolean()) {
						// start with repeat
						f = createFragment(r.begin, fragmentLength, true);
						if (repeats.get(f.getHighBreakend().end) != null) return null; // other end of fragment is also in a repeat
					} else {
						// end with repeat
						f = createFragment(r.end - fragmentLength, fragmentLength, false);
						if (repeats.get(f.getLowBreakend().start) != null) return null; // other end of fragment is also in a repeat
					}
					return f;
				}
			});
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
