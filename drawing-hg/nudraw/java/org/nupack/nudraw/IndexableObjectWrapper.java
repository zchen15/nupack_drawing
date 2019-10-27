package org.nupack.nudraw;

// If an object is going to have an index in multiple collections, rather than have multiple index fields
// in the object, we can use a wrapper in each collection to keep track of the indexes
public class IndexableObjectWrapper<T> extends IndexableObject implements Indexable {

	private T contents;
	
	public IndexableObjectWrapper(T object) {
		this.contents = object;
	}

	public T getContents() {
		return contents;
	}

	public void setContents(T contents) {
		this.contents = contents;
	}

}
