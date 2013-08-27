/**
 * 
 */
package net.wehi.socrates.util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

/**
 * @author hsu
 *
 * Created on Feb 19, 2013
 */
public class GenomicIntervalList implements Iterator<GenomicIndexInterval>, Iterable<GenomicIndexInterval>{
	private ArrayList<GenomicIndexInterval> _backend;
	private boolean finished;
	
	// iterator variables
	private int iterCounter = 0;
	
	public GenomicIntervalList() {
		_backend = new ArrayList<GenomicIndexInterval>();
		finished = false;
	}
	
	public void add(GenomicIndexInterval interval) {
		if (!finished) _backend.add(interval);
	}
	
	public void unique() {
		if (finished) return;
		finished = true;
		Collections.sort(_backend);
		ArrayList<GenomicIndexInterval> uniqueIntervals = new ArrayList<GenomicIndexInterval>();
		for (GenomicIndexInterval interval : _backend) {
			if (uniqueIntervals.size()==0) {
				uniqueIntervals.add(interval);
				continue;
			}
						
			GenomicIndexInterval last = uniqueIntervals.get(uniqueIntervals.size()-1);
			if (last.contains(interval)) continue;
						
			if (last.intersects(interval) || last.adjoins(interval)) {
				last.start = Math.min(last.start, interval.start);
				last.end = Math.max(last.end, interval.end);
				continue;
			}
						
			uniqueIntervals.add( interval );
		}
		_backend = uniqueIntervals;
	}
	
	public int size() { return _backend.size(); }
	
	public boolean hasNext() {
		return (iterCounter < _backend.size());
	}
	
	public GenomicIndexInterval next() {
		return (iterCounter >= _backend.size() ? null : _backend.get(iterCounter++));
	}
	
	public void remove() {
	}
	
	public void reset() {
		iterCounter = 0;
	}
	
	public Iterator<GenomicIndexInterval> iterator() {
		this.reset();
		return this;
	}
}
