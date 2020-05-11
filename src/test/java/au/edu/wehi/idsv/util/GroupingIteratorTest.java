package au.edu.wehi.idsv.util;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertEquals;

public class GroupingIteratorTest {
    @Test
    public void should_group_same_togther() {
        List<List<String>> result = Lists.newArrayList(new GroupingIterator<>(ImmutableList.of("A", "A", "A", "b", "B", "c", "d", "d").iterator(), String.CASE_INSENSITIVE_ORDER));
        assertEquals(4, result.size());
        assertEquals(ImmutableList.of(
                ImmutableList.of("A", "A", "A"),
                ImmutableList.of("b", "B"),
                ImmutableList.of("c"),
                ImmutableList.of("d", "d")), result);
    }
}