package fi.csc.microarray.client.tasks;

import javax.jms.JMSException;

import fi.csc.microarray.databeans.DataManager;

public class LocalTaskExecutor extends TaskExecutor {

	public LocalTaskExecutor(DataManager manager) throws JMSException {
		super(manager);
	}

	@Override
	public void kill(Task task) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public void killAll() {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public void startExecuting(Task task) throws TaskException {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public void startExecuting(Task task, int timeout) throws TaskException {
		throw new UnsupportedOperationException();
	}
}
