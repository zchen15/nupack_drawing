package org.nupack.nudraw;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.log4j.Category;
import org.apache.log4j.Logger;

/**
 * Provides various methods for processing the contents of a Sequence.
 * @author mbp
 *
 */
public class SequenceProcessor {
	
	private static final Category log = Logger.getLogger(ConsoleRunner.class);
	
	private Sequence sequence;
	private List<Helix> helices;
	private List<Loop> loops;
	
	/**
	 * Create a sequence processor for the sequence
	 * @param sequence
	 */
	public SequenceProcessor (Sequence sequence) {
		this.sequence = sequence;
		helices = new AutoIndexingList<Helix>();
		loops = new AutoIndexingList<Loop>();
	}
	
	public void runCalculationChain() {
		linkPairedBases();
		setHelixStructure();
		
		Iterator i = helices.listIterator();
		
		while(i.hasNext()) {
			log.debug(i.next());
		}
	}
	
	/**
	 * Match up corresponding base pairs with a 2-way link
	 * @throws Error
	 */
	private void linkPairedBases() throws Error {
		List<Base> parseStack = new ArrayList<Base>();
		Sequence.BaseIterator bases = sequence.getBaseIterator();
		
		while(bases.hasNext()) {
			Base current = bases.next();
			
			switch (current.getType()) {
			case PAIR_LEFT:
				parseStack.add(current);
				break;
			case PAIR_RIGHT:
				try {
					Base left = parseStack.remove(parseStack.size() - 1);
					left.setPairedBase(current);
					current.setPairedBase(left);
				}
				catch (ArrayIndexOutOfBoundsException e) {
					throw new Error("Parse stack is empty");
				}
				break;
			default:
				break;
			}
		}
		
		if (parseStack.size() > 0) {
			throw new Error("Parse stack is not empty");
		}
	}
	
	/**
	 * Calculates the loop and helix structure
	 *
	 */
	private void setHelixStructure() {
		
		Sequence.BaseIterator bases = sequence.getBaseIterator();
		
		Base b = null;
		
		// There's index - 1 going on in the loop, so rather than have special cases, we'll just run the first by hand
		if (bases.hasNext())
		{
			b = bases.next();
			
			if (b.isPaired()) {
				newHelix().linkBase(b);
			}
		}
		
		while(bases.hasNext()) {
			b = bases.next();
			
			// helices are only applicable to paired bases
			if (!b.isPaired()) {
				continue;
			}
			
			// if paired to something earlier, in the same helix as the paired
			// base
			if (b.getPairedBase().getIndex() < b.getIndex()) {
				b.getPairedBase().getHelix().linkBase(b);
			} else if (sequence.getBase(b.getIndex() - 1).getPairedBase() == null	// prev base unpaired; new helix
					|| sequence.getBase(b.getIndex() - 1).getPairedBase() != 		// prev base paired; bulge loop
						sequence.getBase(b.getPairedBase().getIndex() + 1)
					|| !b.isConnectedPrev()											// strand break on this side of the helix
					|| !b.getPairedBase().isConnectedNext()) {						// strand break on other side
				newHelix().linkBase(b);
			} else {
				sequence.getBase(b.getIndex() - 1).getHelix().linkBase(b); 			// use helix of previous base
			}
		}
	}

	public void setLoopStructure() {
		
	}
	
	public Sequence getSequence() {
		return sequence;
	}

	/* convenience methods */
	
	private Helix newHelix() {
		Helix h = new Helix();
		helices.add(h);
		return h;
	}
	
	private Loop newLoop() {
		Loop l = new Loop();
		loops.add(l);
		return l;
	}

}
