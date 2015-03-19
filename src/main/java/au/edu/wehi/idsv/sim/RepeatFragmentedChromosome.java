package au.edu.wehi.idsv.sim;

import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import au.edu.wehi.idsv.ProcessingContext;

import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeMap;
import com.google.common.collect.TreeRangeSet;

public class RepeatFragmentedChromosome extends FragmentedChromosome {
	private static final Log log = Log.getInstance(RepeatFragmentedChromosome.class);
	public RepeatFragmentedChromosome(ProcessingContext context, String chr, int breakMargin, int seed) {
		super(context, chr, breakMargin, seed);
	}
	private RangeSet<Integer> calcInvalidBreakPositions(int basesBuffer) {
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
	private static RangeMap<Integer, RepeatMaskerRepeat> toLookup(List<RepeatMaskerRepeat> repeatList) throws IOException {
		RangeMap<Integer, RepeatMaskerRepeat> elements = TreeRangeMap.create();
		for (RepeatMaskerRepeat repeat : repeatList) {
			elements.put(Range.closed(repeat.begin, repeat.end), repeat);
		}
		return elements;
	}
	public void assembleSingleSidedRepeats(File fasta, File vcf, boolean includeReference,
			int fragments, int fragmentLength,
			File repeatMaskerOutput, String repeatClassFamily) throws IOException {
		log.info(String.format("Loading %s from %s", repeatClassFamily, repeatMaskerOutput));
		List<RepeatMaskerRepeat> repeatList = RepeatMaskerRepeat.loadRepeatMaskerOutput(repeatMaskerOutput, getChr(), repeatClassFamily);
		RangeMap<Integer, RepeatMaskerRepeat> repeats = toLookup(repeatList);
		RangeSet<Integer> invalid = calcInvalidBreakPositions(margin + fragmentLength);
		Collections.shuffle(repeatList, rng);
		int fragmentsAdded = 0;
		List<Fragment> fragList = new ArrayList<Fragment>();
		for (RepeatMaskerRepeat r : repeatList) {
			// skip repeats we can't use
			if (!invalid.subRangeSet(Range.closedOpen(r.begin, r.end + 1)).isEmpty()) continue;
			// choose a direction
			Fragment f = null;
			if (rng.nextBoolean()) {
				// fragment starts with repeat
				int low = r.begin;
				int high = r.begin + fragmentLength;
				byte[] fragmentSeq = Arrays.copyOfRange(seq, low - 1, high);
				if (repeats.get(high) != null) continue; // other end of fragment is also in a repeat
				f = new Fragment(referenceIndex, low, fragmentSeq, true);
			} else {
				// fragment starts with repeat
				int low = r.end - fragmentLength;
				int high = r.end;
				byte[] fragmentSeq = Arrays.copyOfRange(seq, low - 1, high);
				if (repeats.get(low) != null) continue; // other end of fragment is also in a repeat
				f = new Fragment(referenceIndex, low, fragmentSeq, false);
			}
			fragList.add(f);
			fragmentsAdded++;
			if (fragmentsAdded >= fragments) {
				break;
			}
		}
		log.info(String.format("%d fragments created", fragmentsAdded));
		assemble(fasta, vcf, fragList, includeReference);
	}
}
