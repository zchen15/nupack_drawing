package org.nupack.nudraw;

/**
 * Same as a BaseCollection, except this one automatically indexes
 * bases that are added to it
 * @author mbp
 *
 */
public class AutoIndexingBaseCollection extends BaseCollection {

	/**
	 * In this subclass, we use a different List implementation
	 * to provide the auto-indexing capability
	 *
	 */
	public AutoIndexingBaseCollection() {
		super(new AutoIndexingList<Base>());
	}
	
}
