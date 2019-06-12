package gridss.filter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;

import java.util.List;

/**
 * Filters records only if they are filtered by all filters
 * 
 * @author Daniel Cameron
 *
 */
public class UnionAggregateFilter implements SamRecordFilter {
    private final List<SamRecordFilter> filters;
    public UnionAggregateFilter(final List<SamRecordFilter> filters) {
        this.filters = filters;
    }
    public boolean filterOut(final SAMRecord record) {
        for (final SamRecordFilter filter : filters) {
            if (!filter.filterOut(record)) {
                return false;
            }
        }
        return true;
    }
    public boolean filterOut(final SAMRecord first, final SAMRecord second) {
         for (final SamRecordFilter filter : filters) {
            if (!filter.filterOut(first, second)) {
                return false;
            }
        }
        return true;
    }
}
