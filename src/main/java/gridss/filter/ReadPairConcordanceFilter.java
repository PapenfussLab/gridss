package gridss.filter;

import au.edu.wehi.idsv.ReadPairConcordanceCalculator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;

/**
 * Filters out reads that are not part of either a concordant or discordantly aligned read pair
 * @author Daniel Cameron
 *
 */
public class ReadPairConcordanceFilter implements SamRecordFilter {
	private final ReadPairConcordanceCalculator calc;
	private final boolean includeConcordant;
	private final boolean includeDiscordant;

	public ReadPairConcordanceFilter(ReadPairConcordanceCalculator calc, boolean includeConcordant, boolean includeDiscordant) {
		this.calc = calc;
		this.includeConcordant = includeConcordant;
		this.includeDiscordant = includeDiscordant;
		if (includeConcordant && includeDiscordant) {
			throw new IllegalArgumentException("Filtering both concordant and discordant reads removes all reads");
		}
	}
	public boolean isConcordant(SAMRecord r) {
		return calc.isConcordant(r);
	}
	public boolean isDiscordant(SAMRecord r) {
		return r.getReadPairedFlag() &&
		!r.getReadUnmappedFlag() &&
		!r.getMateUnmappedFlag() &&
		!isConcordant(r);
	}
	private boolean shouldKeep(SAMRecord record) {
		return (includeConcordant && isConcordant(record)) ||
				(includeDiscordant && isDiscordant(record));
	}
	@Override
	public boolean filterOut(SAMRecord record) {
		return !shouldKeep(record);
	}

	@Override
	public boolean filterOut(SAMRecord first, SAMRecord second) {
		return filterOut(first) || filterOut(second);
	}
}
