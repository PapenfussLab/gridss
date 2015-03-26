package au.edu.wehi.idsv.sim;

import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import au.edu.wehi.idsv.ProcessingContext;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;

public class RepeatFragmentedChromosome extends FragmentedChromosome {
	private static final Log log = Log.getInstance(RepeatFragmentedChromosome.class);
	private File repeatMaskerOutput;
	private String repeatClassFamily;
	public RepeatFragmentedChromosome(ProcessingContext context, String chr, int breakMargin, int fragmentLength, File repeatMaskerOutput, String repeatClassFamily, int seed) {
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
						// fragment starts with repeat
						f = createFragment(r.begin, fragmentLength, true);
						if (repeats.get(f.getEndBreakend().end) != null) return null; // other end of fragment is also in a repeat
					} else {
						// fragment ends with repeat
						f = createFragment(r.end - fragmentLength, fragmentLength, true);
						if (repeats.get(f.getStartBreakend().start) != null) return null; // other end of fragment is also in a repeat
					}
					return f;
				}
			});
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
