package fi.csc.microarray.client.dataimport;

/**
 * Abstract class for time consuming import process which is 
 * run on its own thread and which state is informed to user by 
 * using <code>ProgressInformator</code>
 * 
 * @author mkoski
 *
 */
public abstract class RunnableImportProcess implements Runnable {

	/**
	 * Running import process
	 */
	private Thread process;
	
	/**
	 * Progress informator
	 */
	protected ProgressInformator informator;
	
	/**
	 * Creates a new RunnableImportProcess with a specified informator
	 * 
	 * @param informator Informator
	 */
	public RunnableImportProcess(ProgressInformator informator) {
		this.informator = informator;
		informator.setProcess(this);
	}
	
	/**
	 * Passes a message of the state of the process to the informator 
	 * which shows the message to the user in some way.
	 * 
	 * @param message Message of the state of the running process
	 */
	protected void setStateInformationMessage(String message){
		this.informator.setMessage(message);
	}
	
	/**
	 * State value of the running process
	 * 
	 * @param value
	 */
	protected void setStateValue(int value){
		this.informator.setValue(value);
	}
	
	protected ProgressInformator getInformator(){
		return this.informator;
	}
	
	/**
	 * Time consuming task to be done.
	 */
	public abstract void taskToDo();

	public void stopProcess() {
		if(process == null){
			throw new IllegalThreadStateException("The process is not running");
		} else {
			process.interrupt();
		}
	}
	
	/**
	 * Creates an informator and runs an import process on its own thread.
	 * This is only for internal use, programmers should call runProcess() instead.
	 */
	public final void run(){
		this.informator.initializeInformator();
		taskToDo();
		this.informator.destroyInformator();
		this.informator = null;
	}
	
	/**
	 * Run the lines in taskToDo()
	 */
	public void runProcess(){
		process = new Thread(this, "Import thread");
		process.start();
	}
}
