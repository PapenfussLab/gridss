package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.util.List;

import org.junit.Before;
import org.junit.Test;

import com.google.common.base.Function;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;

import htsjdk.samtools.SAMRecord;


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
					getContext().getAssemblyParameters().applyAnnotationFilters(input);
					return input.getSAMRecord();
				} }),
				Iterators.transform(in.iterator(), new Function<SAMRecordAssemblyEvidence, SAMRecord>() {
					@Override
					public SAMRecord apply(SAMRecordAssemblyEvidence input) {
						getContext().getAssemblyParameters().applyAnnotationFilters(input);
						return input.getRemoteSAMRecord();
					} })
			, true, true)));
	}
	@Test
	public void should_filter_reference_breakend() {
		List<DirectedEvidence> evidence = Lists.newArrayList();
		evidence.add(SCE(FWD, Read(0, 5, "5M5S")));
		//evidence.add(SCE(FWD, Read(0, 5, "5M6S")));
		//evidence.add(SCE(FWD, Read(0, 5, "5M7S")));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, Lists.transform(evidence, EID),
				0, 10, 5, B("AAAAAAAAAA"), B("AAAAAAAAAA"));
		e.hydrateEvidenceSet(evidence);
		e.annotateAssembly();
		in.add(e);
		assertEquals(0, out.size());
	}
	@Test
	public void should_filter_reference_breakpoint() {
		List<DirectedEvidence> evidence = Lists.newArrayList();
		evidence.add(SCE(FWD, Read(0, 5, "5M5S")));
		evidence.add(SCE(FWD, Read(0, 5, "5M6S")));
		evidence.add(SCE(FWD, Read(0, 5, "5M7S")));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, Lists.transform(evidence, EID),
				0, 10, 5, B("AAAAAAAAAA"), B("AAAAAAAAAA"));
		e.hydrateEvidenceSet(evidence);
		e.annotateAssembly();
		SAMRecord r = Read(0, 11, "5M");
		r.setReadName(BreakpointFastqEncoding.getRealignmentFastq(e).getReadHeader());
		in.add(AssemblyFactory.incorporateRealignment(getContext(), e, ImmutableList.of(r)));
		go();
		assertEquals(0, out.size());
	}
	@Test
	public void should_ignore_filtered_assembly_breakend() {
		SAMRecordAssemblyEvidence e = new SAMRecordAssemblyEvidenceIteratorTest().BE(1);
		getContext().getAssemblyParameters().applyAnnotationFilters(e);
		assertFalse(new SAMRecordAssemblyEvidenceFilteringIterator(getContext(), Iterators.singletonIterator(e)).hasNext());
	}
	@Test
	public void should_filter_reference_allele() {
		List<DirectedEvidence> evidence = Lists.newArrayList();
		evidence.add(SCE(FWD, Read(0, 10, "5S5M5S")));
		evidence.add(SCE(FWD, Read(0, 10, "6S5M6S")));
		evidence.add(SCE(FWD, Read(0, 10, "7S5M7S")));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), BWD, Lists.transform(evidence, EID),
				0, 10, 5, B("AAAAAAAAAA"), B("AAAAAAAAAA"));
		e.hydrateEvidenceSet(evidence);
		e.annotateAssembly();
		SAMRecord r = Read(0, 5, "5M");
		r.setReadName(BreakpointFastqEncoding.getRealignmentFastq(e).getReadHeader());
		in.add(AssemblyFactory.incorporateRealignment(getContext(), e, ImmutableList.of(r)));
		go();
		assertEquals(0, out.size());
	}
	@Test
	public void should_not_filter_valid_breakpoints() {
		List<DirectedEvidence> evidence = Lists.newArrayList();
		evidence.add(SCE(FWD, Read(0, 10, "5S5M5S")));
		evidence.add(SCE(FWD, Read(0, 10, "6S5M6S")));
		evidence.add(SCE(FWD, Read(0, 10, "7S5M7S")));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), BWD, Lists.transform(evidence, EID),
				0, 10, 5, B("AAAAAAAAAA"), B("AAAAAAAAAA"));
		e.hydrateEvidenceSet(evidence);
		e.annotateAssembly();
		SAMRecord r = Read(0, 105, "5M");
		r.setReadNegativeStrandFlag(true);
		r.setReadName(BreakpointFastqEncoding.getRealignmentFastq(e).getReadHeader());
		in.add(AssemblyFactory.incorporateRealignment(getContext(), e, ImmutableList.of(r)));
		go();
		assertEquals(1, out.size());
	}
}
