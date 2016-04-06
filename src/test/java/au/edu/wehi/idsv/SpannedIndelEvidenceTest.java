package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.SAMRecord;

import java.util.List;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;


public class SpannedIndelEvidenceTest extends TestHelper {	
	@Test
	public void qual_should_use_library_indel_distribution_and_split_score_to_both_sides() {
		SAMRecord r = Read(2, 1, "5M5D5M");
		r.setMappingQuality(40);
		SoftClipEvidenceIterator it = new SoftClipEvidenceIterator(SES(), Iterators.singletonIterator(r));
		SpannedIndelEvidence e = (SpannedIndelEvidence)it.next();
		it.close();
		assertEquals(5, e.getBreakendQual() + e.asRemote().getBreakendQual(), 0.01);
		assertEquals(5, e.getBreakpointQual() + e.asRemote().getBreakpointQual(), 0.01);
	}
	@Test
	public void should_allow_source_read_alignment_to_either_strand() {
		SAMRecord r = Read(2, 1, "2M2D2M");
		r.setMappingQuality(40);
		r.setReadNegativeStrandFlag(true);
		SpannedIndelEvidence e = IE(r);
		assertEquals(new BreakpointSummary(2, FWD, 2, 2, 2, BWD, 5, 5), e.getBreakendSummary());
		assertEquals(new BreakpointSummary(2, BWD, 5, 5, 2, FWD, 2, 2), e.asRemote().getBreakendSummary());
	}
	@Test
	public void should_span_deletion_event() {
		SAMRecord r = Read(2, 1, "2M2D2M");
		r.setReadBases(B("NNNN"));
		r.setMappingQuality(40);
		SpannedIndelEvidence e = IE(r);
		assertEquals(new BreakpointSummary(2, FWD, 2, 2, 2, BWD, 5, 5), e.getBreakendSummary());
		assertEquals(new BreakpointSummary(2, BWD, 5, 5, 2, FWD, 2, 2), e.asRemote().getBreakendSummary());
	}
	@Test
	public void should_span_insertion_event() {
		SAMRecord r = Read(2, 1, "2M2I2M");
		r.setReadBases(B("NNTTNN"));
		r.setMappingQuality(40);
		SpannedIndelEvidence e = IE(r);
		assertEquals(new BreakpointSummary(2, FWD, 2, 2, 2, BWD, 3, 3), e.getBreakendSummary());
		assertEquals(new BreakpointSummary(2, BWD, 3, 3, 2, FWD, 2, 2), e.asRemote().getBreakendSummary());
		assertEquals("TT", e.getUntemplatedSequence());
		assertEquals("TT", e.asRemote().getUntemplatedSequence());
	}
	@Test
	public void indels_should_have_unique_evidenceID() {
		SAMRecord r = Read(2, 1, "2M2I1M2D2M");
		r.setReadBases(B("NNNNNNN"));
		r.setMappingQuality(40);
		List<SpannedIndelEvidence> e = IEList(SES(), r);
		assertEquals(2, e.size());
		assertEquals(4, e.stream().flatMap(ie -> ImmutableList.of(ie, ie.asRemote()).stream()).map(ie -> ie.getEvidenceID()).distinct().count());
	}
}
