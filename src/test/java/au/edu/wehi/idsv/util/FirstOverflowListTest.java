package au.edu.wehi.idsv.util;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.ListIterator;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;


public class FirstOverflowListTest {
	@Test
	public void should_allow_null_alternates() {
		noalt();
	}
	@Test(expected=IllegalArgumentException.class)
	public void should_not_allow_null_primary() {
		new FirstOverflowList<Integer>(null, ImmutableList.of(ImmutableList.of(1)));
	}
	public FirstOverflowList<Integer> noalt() {
		return new FirstOverflowList<Integer>(ImmutableList.of(0, 1), null);
	}
	public FirstOverflowList<Integer> jagged() {
		List<Integer> p = Lists.newArrayList(0, 1, 2, 3);
		List<List<Integer>> a = Lists.newArrayList();
		a.add(ImmutableList.of(4, 5, 6));
		a.add(ImmutableList.of(7, 8));
		a.add(ImmutableList.of(9));
		a.add(null);
		return new FirstOverflowList<Integer>(p, a);
	}
	@Test
	public void child_collections_should_have_primary_as_first_element() {
		assertEquals(Lists.newArrayList(0, 4, 5, 6), jagged().get(0));
		assertEquals(Lists.newArrayList(1, 7, 8   ), jagged().get(1));
		assertEquals(Lists.newArrayList(2, 9      ), jagged().get(2));
		assertEquals(Lists.newArrayList(3         ), jagged().get(3));
		
		assertEquals(Lists.newArrayList(0         ), noalt().get(0));
		assertEquals(Lists.newArrayList(1         ), noalt().get(1));
	}
	@Test
	public void get_2d() {
		assertEquals(0, (int)jagged().get(0, 0));
		assertEquals(4, (int)jagged().get(0, 1));
		assertEquals(5, (int)jagged().get(0, 2));
		assertEquals(6, (int)jagged().get(0, 3));
		assertEquals(1, (int)jagged().get(1, 0));
		assertEquals(7, (int)jagged().get(1, 1));
		assertEquals(8, (int)jagged().get(1, 2));
		assertEquals(2, (int)jagged().get(2, 0));
		assertEquals(9, (int)jagged().get(2, 1));
		assertEquals(3, (int)jagged().get(3, 0));
		
		assertEquals(0, (int)noalt().get(0, 0));
		assertEquals(1, (int)noalt().get(1, 0));
	}
	@Test
	public void size() {
		assertEquals(4, jagged().size());
		assertEquals(2, noalt().size());
	}
	@Test
	public void size_2d() {
		assertEquals(4, jagged().size(0));
		assertEquals(3, jagged().size(1));
		assertEquals(2, jagged().size(2));
		assertEquals(1, jagged().size(3));
		assertEquals(1, noalt().size(0));
		assertEquals(1, noalt().size(1));
	}
	@Test
	public void asFlattenedIterable() {
		assertEquals(Lists.newArrayList(Iterables.concat(jagged())), Lists.newArrayList(jagged().asFlattenedIterable()));
	}
	@Test(expected=IndexOutOfBoundsException.class)
	public void get_1d_should_throw_out_of_bounds() {
		jagged().get(8);
	}
	@Test(expected=IndexOutOfBoundsException.class)
	public void get_2d_should_throw_out_of_bounds_on_null_child() {
		jagged().get(3, 1);
	}
	@Test(expected=IndexOutOfBoundsException.class)
	public void get_2d_should_throw_out_of_bounds() {
		jagged().get(2, 2);
	}
	@Test
	public void ChildListListIterator_prev_next_should_return_same_result() {
		ListIterator<List<Integer>> li = jagged().listIterator(1);
		assertEquals(li.next(), li.previous());
		assertEquals(li.next(), li.previous());
	}
	@Test
	public void ChildListListIterator_nextIndex() {
		ListIterator<List<Integer>> li = jagged().listIterator();
		assertEquals(0, li.nextIndex());
		li.next();
		assertEquals(1, li.nextIndex());
	}
	@Test
	public void ChildListListIterator_previousIndex() {
		ListIterator<List<Integer>> li = jagged().listIterator();
		assertEquals(-1, li.previousIndex());
		li.next();
		li.next();
		assertEquals(1, li.previousIndex());
		li.previous();
		assertEquals(0, li.previousIndex());
	}
}
