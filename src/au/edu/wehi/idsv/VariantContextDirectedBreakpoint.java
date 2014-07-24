package au.edu.wehi.idsv;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.List;

import au.edu.wehi.idsv.util.CollectionUtil;

import com.google.common.base.Function;
import com.google.common.collect.Iterables;
import com.google.common.primitives.Bytes;

public class VariantContextDirectedBreakpoint extends VariantContextDirectedEvidence implements DirectedBreakpoint {
	public VariantContextDirectedBreakpoint(ProcessingContext processContext, EvidenceSource source, VariantContext context) {
		super(processContext, source, context);
		assert(super.getBreakendSummary() instanceof BreakpointSummary);
	}
	@Override
	public BreakpointSummary getBreakendSummary() {
		return (BreakpointSummary)super.getBreakendSummary();
	}
	@Override
	public int getRemoteMapq() {
		return getMapqAssemblyRemoteMax();
	}
	@Override
	public int getRemoteBaseLength() {
		return getAssemblyBreakendLengthMax();
	}
	@Override
	public int getRemoteBaseCount() {
		return getAssemblyBaseCount(null);
	}
	@Override
	public int getRemoteMaxBaseQual() {
		byte[] qual = getBreakendQuality();
		if (qual == null || qual.length == 0) return 0;
		List<Byte> list = Bytes.asList(qual);
		return CollectionUtil.maxInt(Iterables.transform(list, new Function<Byte, Integer>() {
			@Override
			public Integer apply(Byte arg0) {
				return (Integer)(int)(byte)arg0;
			}}), 0);
	}
	@Override
	public int getRemoteTotalBaseQual() {
		// TODO Auto-generated method stub
		return 0;
	}
}
