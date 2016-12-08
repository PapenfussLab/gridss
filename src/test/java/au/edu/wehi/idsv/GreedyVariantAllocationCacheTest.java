package au.edu.wehi.idsv;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.ImmutableList;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;

public class GreedyVariantAllocationCacheTest extends TestHelper {
	@Test
	public void should_return_first_best_variant() {
		DirectedBreakpoint[] evidence = new DirectedBreakpoint[] {
			(DiscordantReadPair)NRRP(DP(0, 1, "1M", true, 1, 1, "1M", true)),
			SR(Read(0, 1, "10M5S"), Read(1, 100, "5M")),
			IE(Read(0, 1, "1M1D1M"))
		};
		GreedyVariantAllocationCache cache = new GreedyVariantAllocationCache(true, true, true);
		for (DirectedBreakpoint e : evidence) {
			cache.addBreakpoint("variant1", 1, e);
			cache.addBreakpoint("variant2", 2, e);
			cache.addBreakpoint("variant3", 2, e);
		}
		for (DirectedBreakpoint e : evidence) {
			assertFalse(cache.isBestBreakpoint("variant1", e));
			assertTrue(cache.isBestBreakpoint("variant2", e));
			assertFalse(cache.isBestBreakpoint("variant3", e));
			assertFalse(cache.isBestBreakpoint("variant4", e));
		}
	}
	@Test
	public void dp_must_be_best_call_for_read_pair() {
		DiscordantReadPair[] evidence = new DiscordantReadPair[] {
				(DiscordantReadPair)NRRP(withReadName("read", DP(1, 1, "1M", true, 2, 2, "1M", true))),
				(DiscordantReadPair)NRRP(withReadName("read", DP(2, 2, "1M", true, 1, 1, "1M", true))),
				(DiscordantReadPair)NRRP(withReadName("read", DP(3, 3, "1M", true, 4, 4, "1M", true))),
				(DiscordantReadPair)NRRP(withReadName("read", DP(1, 1, "1M", true, 4, 4, "1M", true))),
		};
		GreedyVariantAllocationCache cache = new GreedyVariantAllocationCache(true, true, true);
		cache.addBreakpoint("variant1", 2, evidence[0]);
		cache.addBreakpoint("variant1", 1, evidence[1]);
		cache.addBreakpoint("variant2", 1, evidence[2]);
		cache.addBreakpoint("variant3", 1, evidence[3]);
		
		assertTrue(cache.isBestBreakpoint("variant1", evidence[0]));
		assertFalse(cache.isBestBreakpoint("variant2", evidence[0]));
		assertFalse(cache.isBestBreakpoint("variant3", evidence[0]));
		
		assertTrue(cache.isBestBreakpoint("variant1", evidence[1]));
		assertFalse(cache.isBestBreakpoint("variant2", evidence[1]));
		assertFalse(cache.isBestBreakpoint("variant3", evidence[1]));		
		
		assertFalse(cache.isBestBreakpoint("variant1", evidence[2]));
		assertFalse(cache.isBestBreakpoint("variant2", evidence[2]));
		assertFalse(cache.isBestBreakpoint("variant3", evidence[2]));
		// Although it has the same local mapping, the other side is incorrect
		// and this evidence is not part of variant1 
		assertFalse(cache.isBestBreakpoint("variant1", evidence[3]));
		assertFalse(cache.isBestBreakpoint("variant2", evidence[3]));
		assertFalse(cache.isBestBreakpoint("variant3", evidence[3]));
	}
	@Test
	public void sr_must_be_best_split() {
		DirectedBreakpoint[] evidence = new DirectedBreakpoint[] {
				SR(withReadName("read", Read(0, 1, "10M10S"))[0], withReadName("read", Read(0, 10, "10M"))[0]),
				SR(withReadName("read", Read(0, 1, "10M10S"))[0], withReadName("read", Read(0, 20, "10M"))[0]),
				SR(withReadName("read", Read(0, 1, "10S10M"))[0], withReadName("read", Read(0, 10, "10M"))[0]),
		};
		GreedyVariantAllocationCache cache = new GreedyVariantAllocationCache(true, true, true);
		cache.addBreakpoint("variant1", 2, evidence[0]);
		cache.addBreakpoint("variant2", 1, evidence[1]);
		cache.addBreakpoint("variant3", 1, evidence[2]);
		
		assertTrue(cache.isBestBreakpoint("variant1", evidence[0]));
		assertFalse(cache.isBestBreakpoint("variant2", evidence[0]));
		assertFalse(cache.isBestBreakpoint("variant3", evidence[0]));
		
		assertFalse(cache.isBestBreakpoint("variant1", evidence[1]));
		assertFalse(cache.isBestBreakpoint("variant2", evidence[1]));
		assertFalse(cache.isBestBreakpoint("variant3", evidence[1]));		
		
		assertFalse(cache.isBestBreakpoint("variant1", evidence[2]));
		assertFalse(cache.isBestBreakpoint("variant2", evidence[2]));
		assertFalse(cache.isBestBreakpoint("variant3", evidence[2]));
	}
	@Test
	public void indel_must_be_best_alignment() {
		DirectedBreakpoint[] evidence = new DirectedBreakpoint[] {
				IE(SES(), withReadName("read", Read(0, 1, "10M5D10M"))[0], 1),
				IE(SES(), withReadName("read", Read(0, 1, "10M6D10M"))[0], 1),
		};
		GreedyVariantAllocationCache cache = new GreedyVariantAllocationCache(true, true, true);
		cache.addBreakpoint("variant1", 2, evidence[0]);
		cache.addBreakpoint("variant2", 1, evidence[1]);
		
		cache.addBreakpoint("variant1", 2, evidence[0]);
		cache.addBreakpoint("variant2", 1, evidence[1]);
		
		assertTrue(cache.isBestBreakpoint("variant1", evidence[0]));
		assertFalse(cache.isBestBreakpoint("variant2", evidence[0]));
		
		assertFalse(cache.isBestBreakpoint("variant1", evidence[1]));
		assertFalse(cache.isBestBreakpoint("variant2", evidence[1]));
	}
	@Test
	public void indel_must_include_both_sides() {
		IndelEvidence[] evidence = new IndelEvidence[] {
				IE(SES(), withReadName("read", Read(0, 1, "10M5D10M"))[0], 1),
				IE(SES(), withReadName("read", Read(0, 1, "10M6D10M"))[0], 1),
		};
		GreedyVariantAllocationCache cache = new GreedyVariantAllocationCache(true, true, true);
		cache.addBreakpoint("variant1", 2, evidence[0]);
		cache.addBreakpoint("variant2", 1, evidence[1]);
		cache.addBreakpoint("variant1", 2, evidence[0].asRemote());
		cache.addBreakpoint("variant2", 1, evidence[1].asRemote());
		
		assertTrue(cache.isBestBreakpoint("variant1", evidence[0].asRemote()));
		assertFalse(cache.isBestBreakpoint("variant2", evidence[0].asRemote()));
		
		assertFalse(cache.isBestBreakpoint("variant1", evidence[1].asRemote()));
		assertFalse(cache.isBestBreakpoint("variant2", evidence[1].asRemote()));
	}
	@Test
	public void should_allow_multiple_indels_given_same_alignment() {
		IndelEvidence[] evidence = new IndelEvidence[] {
				IE(SES(), withReadName("read", Read(0, 1, "10M5D10M5D10M"))[0], 1),
				IE(SES(), withReadName("read", Read(0, 1, "10M5D10M5D10M"))[0], 2),
		};
		GreedyVariantAllocationCache cache = new GreedyVariantAllocationCache(true, true, true);
		cache.addBreakpoint("variant1", 2, evidence[0]);
		cache.addBreakpoint("variant2", 1, evidence[1]);
		
		assertTrue(cache.isBestBreakpoint("variant1", evidence[0]));
		assertFalse(cache.isBestBreakpoint("variant2", evidence[0]));
		
		// no problem the read supporting both the indels at the same time
		assertFalse(cache.isBestBreakpoint("variant1", evidence[1]));
		assertTrue(cache.isBestBreakpoint("variant2", evidence[1]));
	}
	@Test
	public void should_allow_multiple_split_reads_given_same_alignment() {
		SAMRecord mid = withReadName("read", Read(0, 100, "10S5M8S"))[0];
		List<FastqRecord> ralist = SplitReadIdentificationHelper.getSplitReadRealignments(mid, false);
		SAMRecord pre = Read(0, 50, "10M");
		SAMRecord post = Read(0, 150, "8M");
		pre.setReadName(ralist.get(0).getReadHeader());
		post.setReadName(ralist.get(1).getReadHeader());
		SplitReadIdentificationHelper.convertToSplitRead(mid, ImmutableList.of(pre, post));
		List<SplitReadEvidence> sr = new ArrayList<>();
		sr.addAll(SplitReadEvidence.create(SES(), pre));
		sr.addAll(SplitReadEvidence.create(SES(), mid));
		sr.addAll(SplitReadEvidence.create(SES(), post));
		// and an alternative alignment
		sr.add(SR(withReadName("read", Read(0, 1, "10M10S"))[0], withReadName("read", Read(0, 20, "10M"))[0]));
				
		GreedyVariantAllocationCache cache = new GreedyVariantAllocationCache(true, true, true);
		cache.addBreakpoint("variant1", 1, sr.get(0));
		cache.addBreakpoint("variant1", 1, sr.get(1));
		cache.addBreakpoint("variant2", 2, sr.get(2));
		cache.addBreakpoint("variant2", 1, sr.get(3));
		cache.addBreakpoint("variant3", 1, sr.get(4));
		
		// both split reads are fine
		assertTrue(cache.isBestBreakpoint("variant1", sr.get(0)));
		assertTrue(cache.isBestBreakpoint("variant1", sr.get(1)));
		assertTrue(cache.isBestBreakpoint("variant2", sr.get(2)));
		assertTrue(cache.isBestBreakpoint("variant2", sr.get(3)));
	}
	@Test
	public void should_call_best_sc_alignment() {
		SoftClipEvidence[] evidence = new SoftClipEvidence[] {
				SCE(FWD, withReadName("read", Read(0, 1, "10M5S"))),
				SCE(FWD, withReadName("read", Read(1, 1, "10M5S"))),
		};
		GreedyVariantAllocationCache cache = new GreedyVariantAllocationCache(true, true, true);
		cache.addBreakpoint("variant1", 1, evidence[0]);
		cache.addBreakpoint("variant2", 2, evidence[1]);
		
		assertFalse(cache.isBestBreakpoint("variant1", evidence[0]));
		assertFalse(cache.isBestBreakpoint("variant2", evidence[0]));
		
		assertFalse(cache.isBestBreakpoint("variant1", evidence[1]));
		assertTrue(cache.isBestBreakpoint("variant2", evidence[1]));
	}
	@Test
	public void should_call_best_oea_alignment() {
		NonReferenceReadPair[] evidence = new NonReferenceReadPair[] {
				NRRP(withReadName("read", OEA(0, 1, "1M", true))),
				NRRP(withReadName("read", OEA(0, 2, "1M", true))),
		};
		GreedyVariantAllocationCache cache = new GreedyVariantAllocationCache(true, true, true);
		cache.addBreakpoint("variant1", 1, evidence[0]);
		cache.addBreakpoint("variant2", 2, evidence[1]);
		
		assertFalse(cache.isBestBreakpoint("variant1", evidence[0]));
		assertFalse(cache.isBestBreakpoint("variant2", evidence[0]));
		
		assertFalse(cache.isBestBreakpoint("variant1", evidence[1]));
		assertTrue(cache.isBestBreakpoint("variant2", evidence[1]));
	}
}
