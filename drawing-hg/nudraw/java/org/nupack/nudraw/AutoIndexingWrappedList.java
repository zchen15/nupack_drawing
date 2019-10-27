package org.nupack.nudraw;

import java.util.ArrayList;
import java.util.Collection;

/**
 * Custom AbstractList implementation that sets the index of any inserted element with
 * the position in the array that it's being inserted into
 * @author mbp
 *
 * @param <E>
 */

// Haven't figured out the generics with this class
@Deprecated
public class AutoIndexingWrappedList<E extends Indexable> extends AutoIndexingList<E>{

    private final ArrayList<IndexableObjectWrapper<E>> contents;
    
    AutoIndexingWrappedList() {
    	contents = new ArrayList<IndexableObjectWrapper<E>>();
    }

    @Override
    public E get(int index) {
        return contents.get(index).getContents();
    }

    @Override
    public E set(int index, E element) {
        
        IndexableObjectWrapper<E> wrapper = contents.get(index);
        E oldValue = wrapper.getContents();
        
        wrapper.setContents(element);

        contents.set(index, wrapper);
        
        return oldValue;
    }
    
    public IndexableObjectWrapper<E> getWrapper(int index) {
    	return contents.get(index);
    }

    @Override
    public boolean add(E element) {
    	add(contents.size(), element);
    	
    	// Collection specifies that add should return true
    	return true;
    }
    
    @Override
    public void add(int index, E element) {
        IndexableObjectWrapper<E> wrapper = new IndexableObjectWrapper<E>(element);
    	wrapper.setIndex(index);
    	
    	contents.add(index, wrapper);
    }
    
    @Override
    public boolean addAll(int index, Collection<? extends E> c) {
    	throw new UnsupportedOperationException();
    }
    
}
