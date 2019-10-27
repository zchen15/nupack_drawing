package org.nupack.nudraw;

/**
 * Represents an individual base in a nucleic acid
 * 
 * @author mbp
 */
public class Base extends IndexableObject {

	private BaseType type = null;
	private Base pairedBase = null;

	private boolean connectedNext = false;
	private boolean connectedPrev = false;
	
	private Helix helix = null;
	private Loop loop = null;
	
	public Base(BaseType type, Base pairedBase) {
		this.type = type;
		this.pairedBase = pairedBase;
	}

	public Base() {
	}

	public String toString() {
		return "Base " + getIndex() + " Type " + type
				+ ", paired to " + (pairedBase == null ? "none" : pairedBase.getIndex())
				+ ", helix " + (getHelix() == null ? "none" : getHelix().getIndex());
	}

	public boolean isPaired() {
		return pairedBase != null;
	}
	
	// getters & setters
	public BaseType getType() {
		return type;
	}

	public void setType(BaseType type) {
		this.type = type;
	}

	public Base getPairedBase() {
		return pairedBase;
	}

	public void setPairedBase(Base pairedBase) {
		this.pairedBase = pairedBase;
	}
	
	public boolean isConnectedNext() {
		return connectedNext;
	}

	public void setConnectedNext(boolean connectedNext) {
		this.connectedNext = connectedNext;
	}

	public boolean isConnectedPrev() {
		return connectedPrev;
	}

	public void setConnectedPrev(boolean connectedPrev) {
		this.connectedPrev = connectedPrev;
	}

	public Helix getHelix() {
		return helix;
	}

	public void setHelix(Helix helix) {
		this.helix = helix;
	}

	public Loop getLoop() {
		return loop;
	}

	public void setLoop(Loop loop) {
		this.loop = loop;
	}

}
