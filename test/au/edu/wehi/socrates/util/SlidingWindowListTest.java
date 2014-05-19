package au.edu.wehi.socrates.util;

import static org.junit.Assert.assertEquals;

import org.junit.Before;
import org.junit.Test;

public class SlidingWindowListTest {

	@Before
	public void setUp() throws Exception {
	}
	@Test
	public void getWindowSize_should_return_window_size() {
		assertEquals(2, new SlidingWindowList<Integer>(2).getWindowSize());
		assertEquals(4, new SlidingWindowList<Integer>(4).getWindowSize());
	}
	@Test
	public void size_should_return_size() {
		SlidingWindowList<Integer> list = new SlidingWindowList<Integer>(3);
		assertEquals(0, list.size());
		list.add(0);
		assertEquals(1, list.size());
		list.add(9, 9);
		assertEquals(10, list.size());
	}
	@Test
	public void get_should_return_value_when_in_range() {
		SlidingWindowList<Integer> list = new SlidingWindowList<Integer>(3);
		list.add(0);
		list.add(1);
		list.add(2);
		assertEquals((Integer)0, list.get(0));
		assertEquals((Integer)1, list.get(1));
		assertEquals((Integer)2, list.get(2));
	}
	@Test
	public void get_should_return_null_when_out_of_range() {
		SlidingWindowList<Integer> list = new SlidingWindowList<Integer>(1);
		list.add(0);
		list.add(1);
		list.add(2);
		assertEquals(null, list.get(0));
		assertEquals(null, list.get(1));
	}
	@Test
	public void add_should_invalidate_old_values() {
		SlidingWindowList<Integer> list = new SlidingWindowList<Integer>(16);
		for (int i = 0; i < 16; i++) list.add(i);
		list.add(45, 45);
		for (int i = 0; i < 45; i++) assertEquals(null, list.get(i));
		assertEquals((Integer)45, list.get(45));
	}
	@Test
	public void add_should_invalidate_only_if_required() {
		SlidingWindowList<Integer> list = new SlidingWindowList<Integer>(16);
		for (int i = 0; i < 20; i++) list.add(i);
		list.add(25, 25);
		for (int i = 0; i < 25-16; i++) assertEquals(null, list.get(i));
		for (int i = 25-16+1; i < 20; i++) assertEquals((Integer)i, list.get(i));
		for (int i = 20; i < 25; i++) assertEquals(null, list.get(i));
		assertEquals((Integer)25, list.get(25));
	}
	@Test
	public void set_should_ignore_invalid() {
		SlidingWindowList<Integer> list = new SlidingWindowList<Integer>(1);
		list.add(0);
		list.add(1); 
		assertEquals(null, list.set(0, 0));
		assertEquals(null, list.get(0));
	}
}
