package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import java.util.List;
import java.util.stream.Collectors;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import com.google.common.collect.Lists;

import htsjdk.samtools.SAMRecord;

public abstract class EvidenceIdentifierGeneratorTest extends TestHelper {
	protected EvidenceIdentifierGenerator gen;
	protected ProcessingContext pc;
	protected SAMEvidenceSource ses;
	public abstract EvidenceIdentifierGenerator getGeneratorToTest();
	@Before 
	public void init() {
		gen = getGeneratorToTest();
		pc = getContext();
		pc.setEvidenceIDGenerator(gen);
		ses = SES(pc);
	}
	@Test
	public void should_extract_alignment_unique_from_evidenceid() {
		NonReferenceReadPair rpe = NRRP(ses, withName("readname", DP(0, 1, "5M1D1M4S", true, 1, 1, "10M", false)));
		SAMRecord r = withName("readname", Read(0, 1, "5M1D1M4S"))[0];
		SoftClipEvidence sce = SCE(FWD, ses, r);
		IndelEvidence ie = IE(ses, r);
		assertEquals(gen.extractAlignmentUniqueName(rpe.getEvidenceID()), gen.extractAlignmentUniqueName(sce.getEvidenceID()));
		assertEquals(gen.extractAlignmentUniqueName(rpe.getEvidenceID()), gen.extractAlignmentUniqueName(ie.getEvidenceID()));
		assertNotEquals(gen.extractAlignmentUniqueName(sce.getEvidenceID()), gen.extractAlignmentUniqueName(SCE(FWD, ses, withName("readname", Read(0, 2, "5M1D1M4S"))[0]).getEvidenceID()));
	}
	@Test
	public void should_have_unique_evidenceid() {
		NonReferenceReadPair rpe = NRRP(ses, withName("readname", DP(0, 1, "5M1D1M4S", true, 1, 1, "10M", false)));
		SAMRecord r = withName("readname", Read(0, 1, "5M1D1M4S"))[0];
		SoftClipEvidence sce = SCE(FWD, ses, r);
		IndelEvidence ie = IndelEvidence.create(ses, r, 1);
		SAMRecord r2 = withName("readname", Read(0, 2, "5M1D1M4S"))[0];
		SoftClipEvidence sce2 = SCE(FWD, ses, r2);
		IndelEvidence ie2 = IndelEvidence.create(ses, r2, 1);
		List<String> ids = Lists.newArrayList(
				rpe.getEvidenceID(),
				sce.getEvidenceID(),
				ie.getEvidenceID(),
				ie.asRemote().getEvidenceID(),
				sce2.getEvidenceID(),
				ie2.getEvidenceID(),
				ie2.asRemote().getEvidenceID());
		assertEquals(ids.size(), ids.stream().distinct().count());
	}
	@Test
	public void should_extract_alignment_unique() {
		NonReferenceReadPair rpe = NRRP(ses, withName("readname", DP(0, 1, "5M1D1M4S", true, 1, 1, "10M", false)));
		SAMRecord r = withName("readname", Read(0, 1, "5M1D1M4S"))[0];
		SoftClipEvidence sce = SCE(FWD, ses, r);
		IndelEvidence ie = IndelEvidence.create(ses, r, 1);
		SAMRecord r2 = withName("readname", Read(0, 2, "5M1D1M4S"))[0];
		SoftClipEvidence sce2 = SCE(FWD, ses, r2);
		IndelEvidence ie2 = IndelEvidence.create(ses,  r2, 1);
		List<String> ids = Lists.newArrayList(
				rpe.getEvidenceID(),
				sce.getEvidenceID(),
				ie.getEvidenceID(),
				ie.asRemote().getEvidenceID(),
				sce2.getEvidenceID(),
				ie2.getEvidenceID(),
				ie2.asRemote().getEvidenceID());
		List<String> unique = ids.stream().map(s -> gen.extractAlignmentUniqueName(s)).collect(Collectors.toList());
		assertEquals(2, unique.stream().distinct().count());
	}
	@Test
	public void should_extract_segment_unique() {
		NonReferenceReadPair rpe = NRRP(ses, withName("readname", DP(0, 1, "5M1D1M4S", true, 1, 1, "10M", false)));
		SAMRecord r = withName("readname", Read(0, 1, "5M1D1M4S"))[0];
		SoftClipEvidence sce = SCE(FWD, ses, r);
		IndelEvidence ie = IndelEvidence.create(ses, r, 1);
		SAMRecord r2 = withName("readname", Read(0, 2, "5M1D1M4S"))[0];
		SoftClipEvidence sce2 = SCE(FWD, ses, r2);
		IndelEvidence ie2 = IndelEvidence.create(ses, r2, 1);
		List<String> ids = Lists.newArrayList(
				rpe.getEvidenceID(),
				sce.getEvidenceID(),
				ie.getEvidenceID(),
				ie.asRemote().getEvidenceID(),
				sce2.getEvidenceID(),
				ie2.getEvidenceID(),
				ie2.asRemote().getEvidenceID());
		List<String> unique = ids.stream().map(s -> gen.extractSegmentUniqueName(s)).collect(Collectors.toList());
		assertEquals(1, unique.stream().distinct().count());
	}
	@Test
	public void should_be_unique_for_each_indel() {
		SAMRecord r = withName("readname", Read(0, 1, "5M1D5M1D5M"))[0];
		IndelEvidence ie = IndelEvidence.create(ses, r, 1);
		IndelEvidence ie2 = IndelEvidence.create(ses, r, 3);
		assertNotEquals(ie.getEvidenceID(), ie2.getEvidenceID());
	}
	@Test
	public void extract_alignment_unique_should_match_alignment_unique_from_underlying_read() {
		SAMRecord r = Read(0, 1, "10M10S");
		SoftClipEvidence sc = SCE(FWD, ses, r);
		HashedEvidenceIdentifierGenerator gen = new HashedEvidenceIdentifierGenerator();
		Assert.assertEquals(gen.extractAlignmentUniqueName(gen.getEvidenceID(sc)), gen.getAlignmentUniqueName(r));
	}
	@Test
	public void extract_segement_unique_should_match_segement_unique_from_underlying_read() {
		SAMRecord r = Read(0, 1, "10M10S");
		SoftClipEvidence sc = SCE(FWD, ses, r);
		HashedEvidenceIdentifierGenerator gen = new HashedEvidenceIdentifierGenerator();
		Assert.assertEquals(gen.extractSegmentUniqueName(gen.getEvidenceID(sc)), gen.getSegmentUniqueName(r));
	}
}
