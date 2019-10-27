package org.nupack.nudraw;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import sun.reflect.generics.reflectiveObjects.NotImplementedException;

/**
 * Class for all classes that hold collections of Bases.
 * Confusingly, not actually a Collection (the interface).
 * @author mbp
 *
 */
abstract public class BaseCollection extends IndexableObject {

	/* have to make this protected so that it can be accessed in the subclass's constructor */
	private final List<Base> bases;
	
	public BaseCollection() {
		bases = new ArrayList<Base>();
	}
	
	/**
	 * Set up a 2-way link
	 * @param b
	 */
	public void linkBase(Base b) {
		throw new NotImplementedException();
		//addBase(b);
	}
	
	/**
	 * Pass in the List instance to use
	 * @param list
	 */
	public BaseCollection(List<Base> list) {
		bases = list;
	}
	
	public String toString() {
		String out = this.getClass().toString() + " " + getIndex() + ", bases: ";
		
		Iterator<Base> b = bases.iterator();
		
		while(b.hasNext()) {
			out += b.next().getIndex() + ", ";
		}
		
		return out;
	}
	
	public void addBase(Base b) {
		bases.add(b);
	}
	
	public Base getBase(Integer i) {
		return bases.get(i);
	}
	
	public Integer size() {
		return bases.size();
	}

}
