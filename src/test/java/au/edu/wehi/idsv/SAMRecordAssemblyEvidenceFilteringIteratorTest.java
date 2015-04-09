package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import htsjdk.samtools.SAMRecord;

import java.util.List;
import java.util.Set;

import org.junit.Before;
import org.junit.Test;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;


public class SAMRecordAssemblyEvidenceFilteringIteratorTest extends TestHelper {
	private List<SAMRecordAssemblyEvidence> in;
	private List<SAMRecordAssemblyEvidence> out;
	@Before
	public void setup() {
		in = Lists.newArrayList();
		out = Lists.newArrayList();
	}
	public static ProcessingContext getContext() {
		ProcessingContext pc = TestHelper.getContext();
		return pc;
	}
	public void go() {
		out = Lists.newArrayList(new SAMRecordAssemblyEvidenceFilteringIterator(getContext(),
				new SAMRecordAssemblyEvidenceReadPairIterator(getContext(), AES(),
			Iterators.transform(in.iterator(), new Function<SAMRecordAssemblyEvidence, SAMRecord>() {
				@Override
				public SAMRecord apply(SAMRecordAssemblyEvidence input) {
					return input.getSAMRecord();
				} }),
				Iterators.transform(in.iterator(), new Function<SAMRecordAssemblyEvidence, SAMRecord>() {
					@Override
					public SAMRecord apply(SAMRecordAssemblyEvidence input) {
						return input.getRemoteSAMRecord();
					} })
			, true)));
	}
	@Test
	public void should_filter_reference_breakend() {
		Set<DirectedEvidence> evidence = Sets.<DirectedEvidence>newHashSet();
		evidence.add(SCE(FWD, Read(0, 5, "5M5S")));
		//evidence.add(SCE(FWD, Read(0, 5, "5M6S")));
		//evidence.add(SCE(FWD, Read(0, 5, "5M7S")));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, evidence,
				0, 10, 5, B("AAAAAAAAAA"), B("AAAAAAAAAA"), 5, 6);
		in.add(e);
		assertEquals(0, out.size());
	}
	@Test
	public void should_filter_reference_breakpoint() {
		Set<DirectedEvidence> evidence = Sets.<DirectedEvidence>newHashSet();
		evidence.add(SCE(FWD, Read(0, 5, "5M5S")));
		evidence.add(SCE(FWD, Read(0, 5, "5M6S")));
		evidence.add(SCE(FWD, Read(0, 5, "5M7S")));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, evidence,
				0, 10, 5, B("AAAAAAAAAA"), B("AAAAAAAAAA"), 5, 6);		
		SAMRecord r = Read(0, 11, "5M");
		r.setReadName(BreakpointFastqEncoding.getRealignmentFastq(e).getReadHeader());
		in.add(AssemblyFactory.incorporateRealignment(getContext(), e, r));
		go();
		assertEquals(0, out.size());
	}
	@Test
	public void should_ignore_filtered_assembly_breakend() {
		assertFalse(new SAMRecordAssemblyEvidenceFilteringIterator(getContext(), Iterators.singletonIterator(new SAMRecordAssemblyEvidenceIteratorTest().BE(1))).hasNext());
	}
	@Test
	public void should_filter_reference_allele() {
		Set<DirectedEvidence> evidence = Sets.<DirectedEvidence>newHashSet();
		evidence.add(SCE(FWD, Read(0, 10, "5S5M5S")));
		evidence.add(SCE(FWD, Read(0, 10, "6S5M6S")));
		evidence.add(SCE(FWD, Read(0, 10, "7S5M7S")));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), BWD, evidence,
				0, 10, 5, B("AAAAAAAAAA"), B("AAAAAAAAAA"), 5, 6);		
		SAMRecord r = Read(0, 5, "5M");
		r.setReadName(BreakpointFastqEncoding.getRealignmentFastq(e).getReadHeader());
		in.add(AssemblyFactory.incorporateRealignment(getContext(), e, r));
		go();
		assertEquals(0, out.size());
	}
	@Test
	public void should_not_filter_valid_breakpoints() {
		Set<DirectedEvidence> evidence = Sets.<DirectedEvidence>newHashSet();
		evidence.add(SCE(FWD, Read(0, 10, "5S5M5S")));
		evidence.add(SCE(FWD, Read(0, 10, "6S5M6S")));
		evidence.add(SCE(FWD, Read(0, 10, "7S5M7S")));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), BWD, evidence,
				0, 10, 5, B("AAAAAAAAAA"), B("AAAAAAAAAA"), 5, 6);		
		SAMRecord r = Read(0, 105, "5M");
		r.setReadNegativeStrandFlag(true);
		r.setReadName(BreakpointFastqEncoding.getRealignmentFastq(e).getReadHeader());
		in.add(AssemblyFactory.incorporateRealignment(getContext(), e, r));
		go();
		assertEquals(1, out.size());
	}
}
