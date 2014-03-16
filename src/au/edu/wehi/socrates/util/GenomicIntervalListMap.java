/**
 * 
 */
package au.edu.wehi.socrates.util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

/**
 * @author hsu
 *
 * Created on Feb 19, 2013
 */
public class GenomicIntervalListMap<E> implements Iterator<Values<E>>, Iterable<Values<E>>{
	private ArrayList<Values<E>> _backend;
	private boolean finished;
	
	// iterator variables
	private int iterCounter = 0;
	
	public GenomicIntervalListMap() {
		_backend = new ArrayList<Values<E>>();
		finished = false;
	}
	
	public void add(GenomicIndexInterval interval, E value) {
		if (!finished) _backend.add( new Values<E>(interval, value) );
	}
	
	public void add(GenomicIndexInterval interval, ArrayList<E> value) {
		if (!finished) _backend.add( new Values<E>(interval, value) );
	}
	
	public void unique() {
		if (finished) return;
		finished = true;
		Collections.sort(_backend);
		ArrayList<Values<E>> uniqueIntervals = new ArrayList<Values<E>>();
		for (Values<E> v: _backend) {
			if (uniqueIntervals.size()==0) {
				uniqueIntervals.add(v);
				continue;
			}
						
			Values<E> last = uniqueIntervals.get(uniqueIntervals.size()-1);
			if (last.gi.contains(v.gi)) continue;
						
			if (last.gi.intersects(v.gi) || last.gi.adjoins(v.gi)) {
				last.gi.start = Math.min(last.gi.start, v.gi.start);
				last.gi.end = Math.max(last.gi.end, v.gi.end);
				last.addValues( v.values );
				continue;
			}
						
			uniqueIntervals.add( new Values<E>(v.gi, v.values) );
		}
		_backend = uniqueIntervals;
	}
	
	public int size() { return _backend.size(); }
	
	public boolean hasNext() {
		return (iterCounter < _backend.size());
	}
	
	public Values<E> next() {
		return (iterCounter >= _backend.size() ? null : _backend.get(iterCounter++));
	}
	
	public void remove() {
	}
	
	public void reset() {
		iterCounter = 0;
	}
	
	public Iterator<Values<E>> iterator() {
		this.reset();
		return this;
	}
}

class Values<E> implements Comparable<Values<E>> {
	GenomicIndexInterval gi;
	ArrayList<E> values;

	public Values(int chromIdx, int start, int end) {
		gi = new GenomicIndexInterval(chromIdx, start, end);
		values = new ArrayList<E>();
	}

	public Values(int chromIdx, int start, int end, char strand) {
		gi = new GenomicIndexInterval(chromIdx, start, end, strand);
		values = new ArrayList<E>();
	}
	
	public Values(GenomicIndexInterval g, E value) {
		gi = g;
		values = new ArrayList<E>();
		values.add( value );
	}
	
	public Values(GenomicIndexInterval g, ArrayList<E> value) {
		gi = g;
		values = new ArrayList<E>();
		values.addAll( value );
	}
	
	public void addValues(ArrayList<E> value) {
		values.addAll( value );
	}
	
	public int compareTo(Values<E> other) {
		return gi.compareTo(other.gi);
	}
}
