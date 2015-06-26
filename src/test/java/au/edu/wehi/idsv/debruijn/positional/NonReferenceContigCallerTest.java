package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.Test;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.TestHelper;

import com.google.common.collect.Lists;


public class NonReferenceContigCallerTest extends TestHelper {
	public int maxReadLength(DirectedEvidence... input) {
		return Arrays.stream(input).mapToInt(e -> {
			int length = e.getLocalBaseLength();
			if (e instanceof NonReferenceReadPair) {
				length = Math.max(length, ((NonReferenceReadPair)e).getNonReferenceRead().getReadLength());
			}
			return length;
		}).max().orElse(0);
	}
	private NonReferenceContigCaller caller;
	private EvidenceTracker trackedIt;
	public List<SAMRecordAssemblyEvidence> go(ProcessingContext pc, boolean collapse, DirectedEvidence... input) {
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, Arrays.stream(input).map(e -> (SAMEvidenceSource)e.getEvidenceSource()).collect(Collectors.toList()), null);
		int maxSupportWidth = aes.getMaxConcordantFragmentSize() - aes.getMinConcordantFragmentSize() + 1;
		int maxReadLength = maxReadLength(input);
		int k = pc.getAssemblyParameters().k;
		int maxPathLength = pc.getAssemblyParameters().positionalMaxPathLengthInBases(maxSupportWidth);
		int maxPathCollapseLength = pc.getAssemblyParameters().positionalMaxPathCollapseLengthInBases(maxSupportWidth);
		SupportNodeIterator supportIt = new SupportNodeIterator(k, Arrays.stream(input).iterator(), maxReadLength);
		trackedIt = new EvidenceTracker(supportIt);
		AggregateNodeIterator agIt = new AggregateNodeIterator(trackedIt);
		Iterator<KmerPathNode> pnIt = new PathNodeIterator(agIt, maxSupportWidth, maxPathLength, k);
		if (collapse) {
			pnIt = new PathCollapseIterator(pnIt, k, maxSupportWidth, maxPathLength, maxPathCollapseLength, pc.getAssemblyParameters().maxBaseMismatchForCollapse);
			pnIt = new PathSimplificationIterator(pnIt, maxPathLength, maxSupportWidth);
		}
		caller = new NonReferenceContigCaller(pnIt, 0, maxSupportWidth, maxReadLength, k, aes, trackedIt);
		List<SAMRecordAssemblyEvidence> assemblies = Lists.newArrayList(caller);
		return assemblies;
	}
	@Test
	public void should_call_simple_SC() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 5, "6M6S")));
		List<SAMRecordAssemblyEvidence> output = go(pc, true, sce);
		assertEquals(1, output.size());
		assertEquals("ACGTGGTTAACC", S(output.get(0).getAssemblySequence()));
		assertEquals(new BreakendSummary(0, FWD, 5, 5), output.get(0).getBreakendSummary());
	}
	@Test
	public void should_handle_sc_repeated_kmers() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(FWD, withSequence("TTTTTTTTTTTT", Read(0, 5, "6M6S")));
		List<SAMRecordAssemblyEvidence> output = go(pc, true, sce);
		assertEquals(1, output.size());
		assertEquals("TTTTTTTTTTTT", S(output.get(0).getAssemblySequence()));
	}
}
