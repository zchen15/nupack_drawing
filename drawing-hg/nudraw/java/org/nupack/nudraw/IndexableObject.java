package org.nupack.nudraw;

/** 
 * Parent class of Indexable things (Base, Helix, etc)
 * @author mbp
 *
 */
public class IndexableObject implements Indexable {

	private Integer index = null;
	
	public IndexableObject() {
	}

	public Integer getIndex() {
		return index;
	}

	public void setIndex(Integer index) {
		this.index = index;
	}

}
