package au.edu.wehi.idsv.debruijn;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.stream.Collectors;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;

public class ContigCategorySupportHelper {
	private static final String CIGAR_SEPARATOR = ",";
	/**
	 * Per base contig support per category
	 * @param hasCoverageToBase presence of support
	 * @return
	 */
	public static String asSupportCigars(int contigLength, List<BitSet> hasCoverageToBase) {
		return hasCoverageToBase.stream()
			.map(x -> asSupportCigar(contigLength, x))
			.collect(Collectors.joining(CIGAR_SEPARATOR));
	}
	public static String asSupportCigar(int contigLength, BitSet hasCoverageToBase) {
		List<CigarElement> list = new ArrayList<>();
		int i = 0;
		int runLength = 1;
		boolean value = hasCoverageToBase.get(i);
		while (i < contigLength) {
			while (i + 1 < contigLength && hasCoverageToBase.get(i + 1) == value) {
				i++;
				runLength++;
			}
			list.add(new CigarElement(runLength, value ? CigarOperator.EQ : CigarOperator.X));
			i++;
			value = hasCoverageToBase.get(i);
			runLength = 1;
		}
		return new Cigar(list).toString();
	}
	public static List<BitSet> asPerCategoryReadCoverage(String supportCigars) {
		return Arrays.stream(supportCigars.split(CIGAR_SEPARATOR))
				.map(x -> asReadCoverage(x))
				.collect(Collectors.toList());
	}
	public static BitSet asReadCoverage(String supportCigar) {
		Cigar c = TextCigarCodec.decode(supportCigar);
		BitSet bs = new BitSet(c.getReadLength());
		int i = 0;
		for (CigarElement ce : c.getCigarElements()) {
			bs.set(i, i + ce.getLength(), ce.getOperator() == CigarOperator.EQ);
			i += ce.getLength();
		}
		return bs;
	}
	public static List<Boolean> supportsBreakendBefore(int position, String supportCigars) {
		return asPerCategoryReadCoverage(supportCigars).stream()
				.map(bs -> supportsBreakendBefore(position, bs))
				.collect(Collectors.toList());
	}
	public static List<Boolean> hasSupportInInterval(int start, int end, String supportCigars) {
		return asPerCategoryReadCoverage(supportCigars).stream()
				.map(bs -> hasSupportInInterval(start, end, bs))
				.collect(Collectors.toList());
	}
	public static boolean hasSupportInInterval(int start, int end, BitSet coverage) {
		return coverage != null
				&& coverage.get(start, end + 1).cardinality() > 0;
	}
	public static boolean supportsBreakendBefore(int position, BitSet coverage) {
		return coverage != null
				&& coverage.get(0, position).cardinality() > 0
				&& coverage.get(position, coverage.size()).cardinality() > 0;
	}
}
