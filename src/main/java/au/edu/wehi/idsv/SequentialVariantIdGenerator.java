package au.edu.wehi.idsv;

import java.util.concurrent.atomic.AtomicInteger;

public class SequentialVariantIdGenerator implements VariantIdGenerator {
	private final AtomicInteger offset = new AtomicInteger(0);
	private final String prefix;
	public SequentialVariantIdGenerator(String prefix) {
		this.prefix = prefix;
	}
	public String generate(BreakpointSummary breakpoint) {
		return String.format("%s%d", prefix, offset.incrementAndGet());
	}
}
