package fi.csc.microarray.client.dataimport;

/**
 * An interface to a informator which shows detailed information 
 * of the state of time consuming process. Informator can be an 
 * implementation of dialog or just a console appeneder or whatever.
 * 
 * @author mkoski
 *
 */
public interface ProgressInformator {
	
	/**
	 * Sets message of the state of the process
	 * 
	 * @param message Message text to be shown to user
	 */
	public abstract void setMessage(String message);
	
	/**
	 * Sets state value of the process. Can be for 
	 * exaple percents done of the process or some 
	 * predefined constant of process state.
	 * 
	 * @param state Job state
	 */
	public abstract void setValue(int state);
	
	/**
	 * Set minimum value of the job state.
	 * 
	 * @param min
	 */
	public abstract void setMinimunValue(int min);
	
	/**
	 * Set maximum value of the job state. This could be 100 for 
	 * percentage values or 15000 for line number read from file and so on.
	 * 
	 * @param max
	 */
	public abstract void setMaximumValue(int max);
	
	/**
	 * Initializes the informator just before the job starts.
	 *
	 */
	public abstract void initializeInformator();
	
	/**
	 * Destroys the informator after job has been done
	 *
	 */
	public abstract void destroyInformator();
	
	/**
	 * Stops the currently running process
	 *
	 */
	public void stopProcess();

	/**
	 * Sets process to informators field. This method is called when the 
	 * process object is being created.
	 * 
	 * @param process
	 */
	public abstract void setProcess(RunnableImportProcess process);
	
	public abstract void setIndeterminate(boolean newValue);
	
}
