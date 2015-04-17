package fi.csc.microarray.client.tasks;

import javax.jms.JMSException;

import fi.csc.chipster.tools.ngs.LocalNGSPreprocess;
import fi.csc.microarray.client.Session;
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
	public void startExecuting(final Task task) throws TaskException {
		if (!task.getOperationID().equals("LocalNGSPreprocess.java")) {
			// stupid exception
			throw new UnsupportedOperationException();
		}


		Runnable taskRunnable = new LocalNGSPreprocess(task);
		Session.getSession().getApplication().runBlockingTask("running " + task.getFullName(), taskRunnable);
		
		
	}
	
	@Override
	public void startExecuting(Task task, int timeout) throws TaskException {
		throw new UnsupportedOperationException();
	}
}
