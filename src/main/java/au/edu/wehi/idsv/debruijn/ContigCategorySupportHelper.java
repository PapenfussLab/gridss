package au.edu.wehi.idsv.debruijn;

import com.google.common.primitives.Floats;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ContigCategorySupportHelper {

	public static int[] packInts(List<int[]> unpacked) {
		return unpacked.stream()
				.flatMapToInt(IntStream::of)
				.toArray();
	}
	public static List<int[]> unpackInts(int[] packed, int length) {
		if (packed == null) return new ArrayList<>();
		assert(packed.length % length == 0);
		List<int[]> unpacked = new ArrayList<>();
		for (int offset = 0; offset + length <= packed.length; offset += length) {
			unpacked.add(Arrays.copyOfRange(packed, offset, offset + length));
		}
		return unpacked;
	}
	public static float[] packFloats(List<float[]> unpacked) {
		return Floats.toArray(unpacked.stream()
				.flatMap(l -> Floats.asList(l).stream())
				.collect(Collectors.toList()));
	}
	public static List<float[]> unpackFloats(float[] packed, int length) {
		if (packed == null) return new ArrayList<>();
		assert(packed.length % length == 0);
		List<float[]> unpacked = new ArrayList<>();
		for (int offset = 0; offset + length <= packed.length; offset += length) {
			unpacked.add(Arrays.copyOfRange(packed, offset, offset + length));
		}
		return unpacked;
	}
}
