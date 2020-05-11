package au.edu.wehi.idsv.util;

import com.google.common.collect.Lists;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertEquals;

public class BufferedIteratorTest {
	@Test
	public void should_preserve_iterator_ordering() {
		List<Integer> list = Lists.newArrayList(1, 2, 3, 4);
		List<Integer> result = Lists.newArrayList(new BufferedIterator<>(list.iterator(), 2));
		assertEquals(list, result);
	}
}
