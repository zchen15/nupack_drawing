package org.nupack.nudraw;

public class Loop extends BaseCollection {

	@Override
	public void linkBase(Base b) {
		b.setLoop(this);
		addBase(b);
	}
}
