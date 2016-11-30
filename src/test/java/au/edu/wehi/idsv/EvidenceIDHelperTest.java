package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import java.util.List;
import java.util.stream.Collectors;

import org.junit.Test;

import com.google.common.collect.Lists;

import htsjdk.samtools.SAMRecord;

public class EvidenceIDHelperTest extends TestHelper {
	@Test
	public void should_extract_alignment_unique_from_evidenceid() {
		NonReferenceReadPair rpe = NRRP(withName("readname", DP(0, 1, "5M1D1M4S", true, 1, 1, "10M", false)));
		SAMRecord r = withName("readname", Read(0, 1, "5M1D1M4S"))[0];
		SoftClipEvidence sce = SCE(FWD, r);
		IndelEvidence ie = IE(r);
		assertEquals(EvidenceIDHelper.extractAlignmentUniqueName(rpe.getEvidenceID()), EvidenceIDHelper.extractAlignmentUniqueName(sce.getEvidenceID()));
		assertEquals(EvidenceIDHelper.extractAlignmentUniqueName(rpe.getEvidenceID()), EvidenceIDHelper.extractAlignmentUniqueName(ie.getEvidenceID()));
		assertNotEquals(EvidenceIDHelper.extractAlignmentUniqueName(sce.getEvidenceID()), EvidenceIDHelper.extractAlignmentUniqueName(SCE(FWD, withName("readname", Read(0, 2, "5M1D1M4S"))[0]).getEvidenceID()));
	}
	@Test
	public void should_have_unique_evidenceid() {
		NonReferenceReadPair rpe = NRRP(withName("readname", DP(0, 1, "5M1D1M4S", true, 1, 1, "10M", false)));
		SAMRecord r = withName("readname", Read(0, 1, "5M1D1M4S"))[0];
		SoftClipEvidence sce = SCE(FWD, r);
		IndelEvidence ie = IndelEvidence.create(SES(), r, 1);
		SAMRecord r2 = withName("readname", Read(0, 2, "5M1D1M4S"))[0];
		SoftClipEvidence sce2 = SCE(FWD, r2);
		IndelEvidence ie2 = IndelEvidence.create(SES(), r2, 1);
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
		NonReferenceReadPair rpe = NRRP(withName("readname", DP(0, 1, "5M1D1M4S", true, 1, 1, "10M", false)));
		SAMRecord r = withName("readname", Read(0, 1, "5M1D1M4S"))[0];
		SoftClipEvidence sce = SCE(FWD, r);
		IndelEvidence ie = IndelEvidence.create(SES(), r, 1);
		SAMRecord r2 = withName("readname", Read(0, 2, "5M1D1M4S"))[0];
		SoftClipEvidence sce2 = SCE(FWD, r2);
		IndelEvidence ie2 = IndelEvidence.create(SES(), r2, 1);
		List<String> ids = Lists.newArrayList(
				rpe.getEvidenceID(),
				sce.getEvidenceID(),
				ie.getEvidenceID(),
				ie.asRemote().getEvidenceID(),
				sce2.getEvidenceID(),
				ie2.getEvidenceID(),
				ie2.asRemote().getEvidenceID());
		List<String> unique = ids.stream().map(s -> EvidenceIDHelper.extractAlignmentUniqueName(s)).collect(Collectors.toList());
		assertEquals(2, unique.stream().distinct().count());
	}
	@Test
	public void should_extract_segment_unique() {
		NonReferenceReadPair rpe = NRRP(withName("readname", DP(0, 1, "5M1D1M4S", true, 1, 1, "10M", false)));
		SAMRecord r = withName("readname", Read(0, 1, "5M1D1M4S"))[0];
		SoftClipEvidence sce = SCE(FWD, r);
		IndelEvidence ie = IndelEvidence.create(SES(), r, 1);
		SAMRecord r2 = withName("readname", Read(0, 2, "5M1D1M4S"))[0];
		SoftClipEvidence sce2 = SCE(FWD, r2);
		IndelEvidence ie2 = IndelEvidence.create(SES(), r2, 1);
		List<String> ids = Lists.newArrayList(
				rpe.getEvidenceID(),
				sce.getEvidenceID(),
				ie.getEvidenceID(),
				ie.asRemote().getEvidenceID(),
				sce2.getEvidenceID(),
				ie2.getEvidenceID(),
				ie2.asRemote().getEvidenceID());
		List<String> unique = ids.stream().map(s -> EvidenceIDHelper.extractSegmentUniqueName(s)).collect(Collectors.toList());
		assertEquals(1, unique.stream().distinct().count());
	}
	@Test
	public void should_use_hash_separator() {
		SAMRecord r = withName("readname", Read(0, 1, "5M1D1M4S"))[0];
		SoftClipEvidence sce = SCE(FWD, r);
		assertEquals("readname#0", EvidenceIDHelper.extractSegmentUniqueName(sce.getEvidenceID()));
		assertEquals("readname#0#polyA#1#+#5M1D1M4S", EvidenceIDHelper.extractAlignmentUniqueName(sce.getEvidenceID()));
	}
	@Test
	public void should_be_unique_for_each_indel() {
		SAMRecord r = withName("readname", Read(0, 1, "5M1D5M1D5M"))[0];
		IndelEvidence ie = IndelEvidence.create(SES(), r, 1);
		IndelEvidence ie2 = IndelEvidence.create(SES(), r, 3);
		assertNotEquals(ie.getEvidenceID(), ie2.getEvidenceID());
	}
}
