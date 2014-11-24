package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;

import java.util.List;
import java.util.Set;

import org.junit.Before;
import org.junit.Test;

import au.edu.wehi.idsv.vcf.VcfFilter;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;


public class SAMRecordAssemblyEvidenceIteratorTest extends TestHelper {
	private List<SAMRecord> realigned;	
	private List<SAMRecordAssemblyEvidence> in;
	private List<SAMRecordAssemblyEvidence> out;
	@Before
	public void setup() {
		realigned = Lists.newArrayList();
		in = Lists.newArrayList();
		out = Lists.newArrayList();
	}
	public static ProcessingContext getContext() {
		ProcessingContext pc = TestHelper.getContext();
		pc.getAssemblyParameters().minReads = 0;
		return pc;
	}
	public void go(boolean includeFiltered) { go(getContext(), includeFiltered); }
	public void go(ProcessingContext pc, boolean includeFiltered) {
		out = Lists.newArrayList(new SAMRecordAssemblyEvidenceIterator(pc, AES(),
				Iterators.transform(in.iterator(), new Function<SAMRecordAssemblyEvidence, SAMRecord>() {
					@Override
					public SAMRecord apply(SAMRecordAssemblyEvidence input) {
						return input.getSAMRecord();
					} })
				, realigned.iterator(), includeFiltered));
		// check output is in order
		//for (int i = 0; i < out.size() - 1; i++) {
		//	BreakendSummary l0 = out.get(i).getBreakendSummary();
		//	BreakendSummary l1 = out.get(i).getBreakendSummary();
		//	assertTrue(l0.referenceIndex < l1.referenceIndex || (l0.referenceIndex == l1.referenceIndex && l0.start <= l1.start));
		//}
	}
	public SAMRecordAssemblyEvidence BE(int position) {
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchored(getContext(), AES(), BWD, null,
				0, position, 1 , B("AA"), new byte[] { 7,7 }, 5, 6);
		return e;
	}
	@Test
	public void should_return_assembly() {
		in.add(BE(1));
		go(true);
		assertEquals(1, out.size());
		assertTrue(out.get(0) instanceof SAMRecordAssemblyEvidence);
	}
	@Test
	public void should_match_assembly_with_realign() {
		in.add(BE(1));
		SAMRecord r = Read(1, 10, "1M");
		r.setReadName("0#1#" + in.get(0).getEvidenceID());
		realigned.add(r);
		go(true);
		assertEquals(1, out.size());
		assertTrue(out.get(0) instanceof RealignedSAMRecordAssemblyEvidence);
		assertTrue(out.get(0).getBreakendSummary() instanceof BreakpointSummary);
	}
	@Test
	public void should_flag_assembly_if_realign_unmapped() {
		in.add(BE(1));
		SAMRecord r = Unmapped(1);
		r.setReadName("0#1#" + in.get(0).getEvidenceID());
		realigned.add(r);
		go(true);
		assertEquals(1, out.size());
		assertFalse(out.get(0) instanceof RealignedSAMRecordAssemblyEvidence);
		assertFalse(out.get(0).getBreakendSummary() instanceof BreakpointSummary);
	}
	@Test
	public void should_allow_realign_in_order_at_same_position() {
		SAMRecord f = withReadName("0#10#fReadName", Read(0, 1, "5M"))[0];
		SAMRecord b = withReadName("0#1#bReadName", Read(0, 1, "5M"))[0];
		in.add(BE(1));
		SAMRecord assemblyRealigned = withReadName("0#1#" + in.get(0).getEvidenceID(), Read(1, 10, "1M"))[0];
		realigned.add(b);
		realigned.add(assemblyRealigned);
		realigned.add(f);
		go(true);
		assertEquals(1, out.size());
		assertTrue(out.get(0) instanceof RealignedSAMRecordAssemblyEvidence);
		assertTrue(out.get(0).getBreakendSummary() instanceof BreakpointSummary);
	}
	@Test
	public void should_require_realign_in_call_position_order() {
		in.add(BE(2));
		realigned.add(withReadName("0#1#bReadName", Read(0, 1, "5M"))[0]);
		realigned.add(withReadName("0#2#" + in.get(0).getEvidenceID(), Read(1, 10, "1M"))[0]);
		realigned.add(withReadName("0#10#fReadName", Read(0, 1, "5M"))[0]);
		go(true);
		assertEquals(1, out.size());
		assertTrue(out.get(0) instanceof RealignedSAMRecordAssemblyEvidence);
		assertTrue(out.get(0).getBreakendSummary() instanceof BreakpointSummary);
	}
	@Test
	public void should_ignore_filtered_assembly_breakend() {
		in.add(BE(1));
		in.get(0).filterAssembly(VcfFilter.ASSEMBLY_REF);
		go(false);
		assertEquals(0, out.size());
	}
	@Test
	public void should_filter_reference_breakpoint() {
		ProcessingContext pc = getContext();
		Set<DirectedEvidence> evidence = Sets.<DirectedEvidence>newHashSet();
		evidence.add(SCE(FWD, Read(0, 5, "5M5S")));
		evidence.add(SCE(FWD, Read(0, 5, "5M6S")));
		evidence.add(SCE(FWD, Read(0, 5, "5M7S")));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchored(pc, AES(), FWD, evidence,
				0, 10, 5, B("AAAAAAAAAA"), B("AAAAAAAAAA"), 5, 6);		
		SAMRecord r = Read(1, 11, "5M");
		r.setReadName(BreakpointFastqEncoding.getRealignmentFastq(e).getReadHeader());
		in.add(e);
		go(pc, false);
		assertEquals("precondition - breakend should not be filtered", 1, out.size());
		setup(); // reset
		in.add(e);
		realigned.add(r);
		go(pc, false);
		assertEquals(0, out.size());
	}
}
