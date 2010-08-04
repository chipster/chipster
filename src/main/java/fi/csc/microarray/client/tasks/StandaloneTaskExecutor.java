package fi.csc.microarray.client.tasks;

import java.beans.PropertyChangeListener;
import java.util.Collection;

import javax.jms.JMSException;

import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.databeans.DataManager;

public class StandaloneTaskExecutor extends TaskExecutor {

	protected StandaloneTaskExecutor(DataManager manager) throws JMSException {
		super(manager);
	}

	@Override
	public void addChangeListener(PropertyChangeListener listener) {
		// TODO Auto-generated method stub
		super.addChangeListener(listener);
	}
	
	@Override
	protected void addToRunningTasks(Task task) {
		// TODO Auto-generated method stub
		super.addToRunningTasks(task);
	}
	
	@Override
	public Task createTask(Operation operation) {
		// TODO Auto-generated method stub
		return super.createTask(operation);
	}
	
	@Override
	public void execute(Task task) throws TaskException {
		// TODO Auto-generated method stub
		super.execute(task);
	}
	
	@Override
	public int getRunningTaskCount() {
		// TODO Auto-generated method stub
		return super.getRunningTaskCount();
	}
	
	@Override
	public Collection<Task> getTasks(boolean onlyRunning, boolean showHidden) {
		// TODO Auto-generated method stub
		return super.getTasks(onlyRunning, showHidden);
	}
	
	@Override
	public boolean isEventsEnabled() {
		// TODO Auto-generated method stub
		return super.isEventsEnabled();
	}
	
	@Override
	public void kill(Task task) {
		// TODO Auto-generated method stub
		super.kill(task);
	}
	
	@Override
	public void killAll() {
		// TODO Auto-generated method stub
		super.killAll();
	}
	
	@Override
	protected void removeFromRunningTasks(Task task) {
		// TODO Auto-generated method stub
		super.removeFromRunningTasks(task);
	}
	
	@Override
	public void setEventsEnabled(boolean eventsEnabled) {
		// TODO Auto-generated method stub
		super.setEventsEnabled(eventsEnabled);
	}
	
	@Override
	public void startExecuting(Task task) throws TaskException {
		// TODO Auto-generated method stub
		super.startExecuting(task);
	}
	
	@Override
	public void startExecuting(Task task, int timeout) throws TaskException {
		// TODO Auto-generated method stub
		super.startExecuting(task, timeout);
	}
	
	

}
