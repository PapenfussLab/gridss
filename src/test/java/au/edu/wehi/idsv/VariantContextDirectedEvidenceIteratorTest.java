package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.util.ArrayList;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

import com.google.common.collect.Lists;


public class VariantContextDirectedEvidenceIteratorTest extends TestHelper {
	private List<SAMRecord> realigned;	
	private List<VariantContext> vcf;
	private List<VariantContextDirectedEvidence> out;
	@Before
	public void setup() {
		realigned = new ArrayList<>();
		vcf = new ArrayList<>();
		out = new ArrayList<>();
	}
	public static ProcessingContext getContext() {
		ProcessingContext pc = TestHelper.getContext();
		pc.getAssemblyParameters().minReads = 0;
		return pc;
	}
	public void go() {
		
		out = Lists.newArrayList(new VariantContextDirectedEvidenceIterator(getContext(), AES(), vcf.iterator(), realigned.iterator()));
		// check output is in order
		//for (int i = 0; i < out.size() - 1; i++) {
		//	BreakendSummary l0 = out.get(i).getBreakendSummary();
		//	BreakendSummary l1 = out.get(i).getBreakendSummary();
		//	assertTrue(l0.referenceIndex < l1.referenceIndex || (l0.referenceIndex == l1.referenceIndex && l0.start <= l1.start));
		//}
	}
	public VariantContextDirectedEvidence BE(int position) {
		return new AssemblyBuilder(getContext(), AES())
			.assemblerName("test")
			.assemblyBases(B("AA"))
			.anchorLength(1)
			.direction(BWD)
			.referenceAnchor(0, position)
			.assembledBaseCount(5, 6)
			.assemblyBaseQuality(new byte[] { 7,7 } )
			.makeVariant();
	}
	@Test
	public void should_return_assembly() {
		VariantContextDirectedEvidence assembly = BE(1);
		vcf.add(new VariantContextBuilder(assembly).make());
		go();
		assertEquals(1, out.size());
		assertTrue(out.get(0) instanceof VariantContextDirectedEvidence);
	}
	@Test
	public void should_match_assembly_with_realign() {
		VariantContextDirectedEvidence assembly = BE(1);
		vcf.add(new VariantContextBuilder(assembly).make());
		SAMRecord r = Read(1, 10, "1M");
		r.setReadName("0#1#" + vcf.get(0).getID());
		realigned.add(r);
		go();
		assertEquals(1, out.size());
		assertTrue(out.get(0) instanceof VariantContextDirectedEvidence);
		assertTrue(out.get(0).getBreakendSummary() instanceof BreakpointSummary);
	}
	@Test
	public void should_flag_assembly_if_realign_unmapped() {
		VariantContextDirectedEvidence assembly = BE(1);
		vcf.add(new VariantContextBuilder(assembly).make());
		SAMRecord r = Unmapped(1);
		r.setReadName("0#1#" + vcf.get(0).getID());
		realigned.add(r);
		go();
		assertEquals(1, out.size());
		assertTrue(out.get(0) instanceof VariantContextDirectedEvidence);
		assertTrue(out.get(0).getBreakendSummary() instanceof BreakendSummary);
	}
	@Test
	public void should_allow_realign_in_order_at_same_position() {
		SAMRecord f = withReadName("0#10#fReadName", Read(0, 1, "5M"))[0];
		SAMRecord b = withReadName("0#1#bReadName", Read(0, 1, "5M"))[0];
		VariantContextDirectedEvidence assembly = BE(1);
		vcf.add(new VariantContextBuilder(assembly).make());
		SAMRecord assemblyRealigned = withReadName("0#1#" + vcf.get(0).getID(), Read(1, 10, "1M"))[0];
		realigned.add(b);
		realigned.add(assemblyRealigned);
		realigned.add(f);
		go();
		assertEquals(1, out.size());
		assertTrue(out.get(0) instanceof VariantContextDirectedBreakpoint);
		assertTrue(out.get(0).getBreakendSummary() instanceof BreakpointSummary);
	}
	@Test
	public void should_require_realign_in_call_position_order() {
		VariantContextDirectedEvidence assembly = BE(2);
		vcf.add(new VariantContextBuilder(assembly).make());
		realigned.add(withReadName("0#1#bReadName", Read(0, 1, "5M"))[0]);
		realigned.add(withReadName("0#2#" + vcf.get(0).getID(), Read(1, 10, "1M"))[0]);
		realigned.add(withReadName("0#10#fReadName", Read(0, 1, "5M"))[0]);
		go();
		assertEquals(1, out.size());
		assertTrue(out.get(0) instanceof VariantContextDirectedBreakpoint);
		assertTrue(out.get(0).getBreakendSummary() instanceof BreakpointSummary);
	}
	@Test
	public void should_ignore_filtered_variants() {
		VariantContextDirectedEvidence assembly = (VariantContextDirectedEvidence)minimalBreakend()
				.filter("FILTERED")
				.make();
		vcf.add(assembly);
		go();
		assertEquals(0, out.size());
	}
}
