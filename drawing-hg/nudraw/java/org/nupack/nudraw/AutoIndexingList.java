package org.nupack.nudraw;

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collection;

/**
 * Custom AbstractList implementation that sets the index of any inserted element with
 * the position in the array that it's being inserted into
 * @author mbp
 *
 * @param <E>
 */
public class AutoIndexingList<E extends Indexable> extends AbstractList<E> {

    private final ArrayList<E> contents;
    
    AutoIndexingList() {
    	contents = new ArrayList<E>();
    }

    @Override
    public E get(int index) {
        return contents.get(index);
    }
    
    @Override
    public E set(int index, E element) {
        E oldValue = contents.get(index);
        
        element.setIndex(index);
        
        contents.set(index, element);
        
        return oldValue;
    }

    @Override
    public int size() {
        return contents.size();
    }
    
    @Override
    public boolean add(E e) {
    	e.setIndex(contents.size());
    	
    	return contents.add(e);
    }
    
    @Override
    public void add(int index, E e) {
    	e.setIndex(index);
    	
    	contents.add(index, e);
    }
    
    @Override
    public boolean addAll(int index, Collection<? extends E> c) {
    	throw new UnsupportedOperationException();
    }
    
}
