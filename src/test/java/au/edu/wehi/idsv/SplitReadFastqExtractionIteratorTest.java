package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;

import htsjdk.samtools.SAMRecord;

public class SplitReadFastqExtractionIteratorTest extends TestHelper {
	@Test
	public void should_extract_above_min_soft_clip_length() {
		assertEquals(0, Iterators.size(new SplitReadFastqExtractionIterator(ImmutableList.of(Read(0, 1, "1S10M")).iterator(), false, 2, 0, true, false, false, false, getContext().getEvidenceIDGenerator())));
		assertEquals(1, Iterators.size(new SplitReadFastqExtractionIterator(ImmutableList.of(Read(0, 1, "2S10M")).iterator(), false, 2, 0, true, false, false, false, getContext().getEvidenceIDGenerator())));
	}
	@Test
	public void should_filter_secondary() {
		SAMRecord r = Read(0, 1, "10S10M");
		r.setSecondaryAlignment(true);
		assertEquals(0, Iterators.size(new SplitReadFastqExtractionIterator(ImmutableList.of(r).iterator(), false, 1, 0, false, false, false, false, getContext().getEvidenceIDGenerator())));
		assertEquals(1, Iterators.size(new SplitReadFastqExtractionIterator(ImmutableList.of(r).iterator(), false, 1, 0, true, false, false, false, getContext().getEvidenceIDGenerator())));
	}
	@Test
	public void should_filter_supplemenary() {
		SAMRecord r = Read(0, 1, "10S10M");
		r.setSupplementaryAlignmentFlag(true);
		assertEquals(0, Iterators.size(new SplitReadFastqExtractionIterator(ImmutableList.of(r).iterator(), false, 1, 0, true, false, false, false, getContext().getEvidenceIDGenerator())));
	}
	@Test
	public void should_filter_unmapped() {
		SAMRecord r = Read(0, 1, "10S10M");
		r.setReadUnmappedFlag(true);
		assertEquals(0, Iterators.size(new SplitReadFastqExtractionIterator(ImmutableList.of(r).iterator(), false, 1, 0, true, false, false, false, getContext().getEvidenceIDGenerator())));
	}
	@Test
	public void should_filter_existing_split_reads_with_SA_tag() {
		SAMRecord r = Read(0, 1, "10S10M");
		r.setAttribute("SA", "polyA,100,+,10M10S,0,0");
		assertEquals(0, Iterators.size(new SplitReadFastqExtractionIterator(ImmutableList.of(r).iterator(), false, 1, 0, true, false, false, false, getContext().getEvidenceIDGenerator())));
	}
	@Test
	public void should_not_filter_existing_split_reads_with_SA_tag_if_realignExistingSplitReads_is_true() {
		SAMRecord r = Read(0, 1, "10S10M");
		r.setAttribute("SA", "polyA,100,+,10M10S,0,0");
		assertEquals(1, Iterators.size(new SplitReadFastqExtractionIterator(ImmutableList.of(r).iterator(), false, 1, 0, true, true, false, false, getContext().getEvidenceIDGenerator())));
	}
	@Test
	public void should_filter_under_average_mapq() {
		SAMRecord r = Read(0, 1, "2S1M");
		r.setBaseQualities(new byte[] { 0, 10, 40 } );
		assertEquals(1, Iterators.size(new SplitReadFastqExtractionIterator(ImmutableList.of(r).iterator(), false, 1, 4, true, false, false, false, getContext().getEvidenceIDGenerator())));
		assertEquals(1, Iterators.size(new SplitReadFastqExtractionIterator(ImmutableList.of(r).iterator(), false, 1, 5, true, false, false, false, getContext().getEvidenceIDGenerator())));
		assertEquals(0, Iterators.size(new SplitReadFastqExtractionIterator(ImmutableList.of(r).iterator(), false, 1, 6, true, false, false, false, getContext().getEvidenceIDGenerator())));
	}
}
