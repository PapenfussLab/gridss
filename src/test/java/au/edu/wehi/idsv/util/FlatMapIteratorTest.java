package au.edu.wehi.idsv.util;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import org.junit.Assert;
import org.junit.Test;

import java.util.List;

public class FlatMapIteratorTest {
    @Test
    public void round_trip_with_batching_iterator() {
        for (int i = 1; i < 8; i++) {
            List<Integer> input = ImmutableList.of(1, 2, 3, 4, 5);
            List<Integer> output = Lists.newArrayList(new FlatMapIterator<>(new BatchingIterator<>(input.iterator(), i)));
            Assert.assertEquals(input, output);
        }
    }
}