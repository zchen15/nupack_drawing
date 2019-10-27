package org.nupack.nudraw;

import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Holds multiple strands
 * 
 * @author mbp
 *
 */
public class Sequence {

	private final List<Strand> strands;
	private Integer numBases;
	
	public Sequence() {
		strands = new AutoIndexingList<Strand>();
		strands.add(new Strand());
		
		numBases = 0;
	}
	
	/**
	 * Switch to a new, empty strand
	 *
	 */
	public void newStrand() {
		strands.add(new Strand());
	}
	
	public void addBase(Base base) {
		strands.get(strands.size() - 1).addBase(base);
		
		numBases++;
	}
	
	public Strand getStrand(int i) {
		return strands.get(i);
	}

	/**
	 * 
	 * @param i
	 * @return The base at the index when viewing the sequence as one big collection of bases
	 */
	public Base getBase(int i) {
		if (i > getNumBases() - 1) {
			throw new ArrayIndexOutOfBoundsException(i);
		}
		
		int baseIndex = i;
		int strandIndex = 0;
		while(baseIndex > strands.get(strandIndex).size() - 1) {
			baseIndex -= strands.get(strandIndex).size();
			strandIndex++;
		}
		
		return strands.get(strandIndex).getBase(baseIndex);
	}
	
	/**
	 * 
	 * @return Number of strands in the sequence
	 */
	public int getNumStrands() {
		return strands.size();
	}
	
	/**
	 * 
	 * @return Total number of bases in all strands
	 */
	public int getNumBases() {
		return numBases;
	}
	
	public BaseIterator getBaseIterator() {
		return new BaseIterator();
	}
	
	/**
	 * Special iterator that hides the underlying strands of a sequence
	 * @author mbp
	 *
	 */
	class BaseIterator implements Iterator {
		
		private int strandIndex;
		private int baseIndex;
		
		private BaseIterator() {
			strandIndex = 0;
			baseIndex = -1; // increment immediately before we lookup
		}

		public boolean hasNext() {
			if (! atLastBaseInStrand()) {
				return true;
			}
			
			else {
				// we're at the end of the last strand
				if (atLastStrand()) {
					return false;
				}
				// there's one more strand with something in it
				else if (strands.get(strandIndex + 1).size() > 0) {
					return true;
				}
				// there is one more strand, but it's empty
				else
				{
					return false;
				}
			}
			
		}

		public Base next() {
			// room left in the current strand
			if (! atLastBaseInStrand()) {
				baseIndex++;
				return getCurrentBase();
			}
			// we've used up the current strand, and there is another strand
			else if (! atLastStrand())
			{

				strandIndex++;
				// if the new strand has anything in it, we're good
				if(getCurrentStrand().size() > 0) {
					baseIndex = 0; // effectively a baseIndex++
					return getCurrentBase();
				} 
				// there is a new strand, but it's empty
				else {
					throw new NoSuchElementException();
				}
			}
			// there is not another strand
			else {
				throw new NoSuchElementException();
			}
		}

		public void remove() {
			throw new UnsupportedOperationException();
		}
		
		private Strand getCurrentStrand() {
			return strands.get(strandIndex);
		}
		
		private Base getCurrentBase() {
			return getCurrentStrand().getBase(baseIndex);
		}
		
		private boolean atLastBaseInStrand() {
			return getCurrentStrand().size() - 1 <= baseIndex ? true : false;
		}
		
		private boolean atLastStrand() {
			return strands.size() - 1 <= strandIndex ? true : false;
		}
	}
}
