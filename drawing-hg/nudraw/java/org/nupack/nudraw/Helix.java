package org.nupack.nudraw;

public class Helix extends BaseCollection {

	private Loop startLoop;
	private Loop endLoop;
	
	@Override
	public void linkBase(Base b) {
		b.setHelix(this);
		addBase(b);
	}
	
	public Loop getEndLoop() {
		return endLoop;
	}
	public void setEndLoop(Loop endLoop) {
		this.endLoop = endLoop;
	}
	public Loop getStartLoop() {
		return startLoop;
	}
	public void setStartLoop(Loop startLoop) {
		this.startLoop = startLoop;
	}
	
}
