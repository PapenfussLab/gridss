package au.edu.wehi.idsv;

import static org.junit.Assert.*;

import org.junit.Test;

import com.google.common.collect.ImmutableList;

import htsjdk.samtools.SAMRecord;

public class DirectedEvidenceIteratorTest extends TestHelper {
	private SAMRecord sampleSplitRead() {
		SAMRecord r = Read(2, 1, "7M3S");
		// MMMMMMMSSS
		// SSMMMMMMMM
		//   ^---^ homology
		r.setReadBases(    B("ACCTAGAGGG"));
		r.setBaseQualities(B("1234567890"));
		r.setMappingQuality(40);
		r.setAttribute("SA", "polyA,100,+,2S8M,0,0");
		return r;
	}
	@Test
	public void should_not_return_both_soft_clip_and_split_read_for_same_read() {
		SAMRecord r = sampleSplitRead();
		ImmutableList<DirectedEvidence> list = ImmutableList.copyOf(new DirectedEvidenceIterator(ImmutableList.of(r).iterator(), null));
		assertEquals(1, list.size());
	}
	@Test
	public void should_ignore_unmapped_reads() {
		SAMRecord r = sampleSplitRead();
		r.setReadUnmappedFlag(true);
		ImmutableList<DirectedEvidence> list = ImmutableList.copyOf(new DirectedEvidenceIterator(ImmutableList.of(r).iterator(), null));
		assertEquals(0, list.size());
	}
	@Test
	public void should_return_OEA() {
		SAMRecord r = OEA(0,  1, "10M", true)[0];
		ImmutableList<DirectedEvidence> list = ImmutableList.copyOf(new DirectedEvidenceIterator(ImmutableList.of(r).iterator(), SES()));
		assertEquals(1, list.size());
		assertEquals(UnmappedMateReadPair.class, list.get(0).getClass());
	}
	@Test
	public void should_return_DP() {
		SAMRecord r = DP(0, 1, "10M", true, 1, 10, "10M", false)[0];
		ImmutableList<DirectedEvidence> list = ImmutableList.copyOf(new DirectedEvidenceIterator(ImmutableList.of(r).iterator(), SES()));
		assertEquals(1, list.size());
		assertEquals(DiscordantReadPair.class, list.get(0).getClass());
	}
	@Test
	public void should_return_indel() {
		SAMRecord r = Read(0, 1, "10M5I5D10M5I10M");
		ImmutableList<DirectedEvidence> list = ImmutableList.copyOf(new DirectedEvidenceIterator(ImmutableList.of(r).iterator(), null));
		assertEquals(4, list.size());
		assertEquals(IndelEvidence.class, list.get(0).getClass());
	}
}
