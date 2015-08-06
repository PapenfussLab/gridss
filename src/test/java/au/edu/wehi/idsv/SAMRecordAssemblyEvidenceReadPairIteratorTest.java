package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

import com.google.common.base.Function;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;


public class SAMRecordAssemblyEvidenceReadPairIteratorTest extends TestHelper {	
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
	public void go(boolean includeFiltered) { go(getContext(), includeFiltered); }
	public void go(ProcessingContext pc, boolean includeFiltered) {
		out = Lists.newArrayList(new SAMRecordAssemblyEvidenceReadPairIterator(pc, AES(),
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
			, true));
	}
	@Test
	public void should_construct_assembly() {
		List<DirectedEvidence> evidence = Lists.newArrayList();
		evidence.add(SCE(FWD, Read(0, 5, "5M5S")));
		evidence.add(SCE(FWD, Read(0, 5, "5M6S")));
		evidence.add(SCE(FWD, Read(0, 5, "5M7S")));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, Lists.transform(evidence, EID),
				0, 10, 5, B("AAAAAAAAAA"), B("AAAAAAAAAA"), new int[] {5, 6});
		in.add(e);
		go(true);
		assertEquals(1, out.size());
	}
	@Test
	public void shoould_load_remote_evidence() {
		SAMRecordAssemblyEvidence e = AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(0, BWD, 1, 10), null, B("AA"), B("AA"), new int[] {2, 0});
		RealignedSAMRecordAssemblyEvidence re = (RealignedSAMRecordAssemblyEvidence)AssemblyFactory.incorporateRealignment(getContext(), e, ImmutableList.of(
				withReadName("0#1#0#ReadName", Read(1, 100, "2M"))[0]));
		assertEquals(new BreakpointSummary(0, BWD, 1, 10, 1, FWD, 101, 110), re.getBreakendSummary());
		assertEquals(new BreakpointSummary(1, FWD, 101, 110, 0, BWD, 1, 10), re.asRemote().getBreakendSummary());
		
		ArrayList<SAMRecord> inlist = new ArrayList<SAMRecord>();
		inlist.add(re.getSAMRecord());
		inlist.add(re.getRemoteSAMRecord());
		ArrayList<SAMRecord> matelist = new ArrayList<SAMRecord>(inlist);
		Collections.reverse(matelist);
		SAMRecordAssemblyEvidenceReadPairIterator it = new SAMRecordAssemblyEvidenceReadPairIterator(getContext(), AES(), inlist.iterator(), matelist.iterator(), true);
		out = Lists.newArrayList(it);
		
		assertEquals(new BreakpointSummary(0, BWD, 1, 10, 1, FWD, 101, 110), out.get(0).getBreakendSummary());
		assertEquals(new BreakpointSummary(1, FWD, 101, 110, 0, BWD, 1, 10), out.get(1).getBreakendSummary());
	}
}
