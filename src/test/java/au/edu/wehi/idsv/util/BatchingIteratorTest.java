package au.edu.wehi.idsv.util;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

public class BatchingIteratorTest {

    @Test
    public void should_group_by_size() {
        ArrayList<List<Integer>> result = Lists.newArrayList(new BatchingIterator<Integer>(ImmutableList.of(1, 2, 3, 4, 5).iterator(), 2));
        Assert.assertEquals(3, result.size());
        Assert.assertEquals(2, result.get(0).size());
        Assert.assertEquals(2, result.get(1).size());
        Assert.assertEquals(1, result.get(2).size());
    }
}