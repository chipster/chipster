package fi.csc.microarray.comp;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.TimeUnit;

import org.apache.log4j.Logger;
/**
 * A process pool for pooling for example R processes.
 * 
 *  
 * The maximum number of processes controlled by this pool may temporarily be exceeded, 
 * due to the synchronization mechanisms used.
 * 
 *
 */
public class ProcessPool {

	private static int WAIT_TIME_BEFORE_FORCE_DESTROY = 10;
	private static TimeUnit WAIT_TIME_BEFORE_FORCE_DESTROY_UNIT = TimeUnit.SECONDS;
	
	static final Logger logger = Logger.getLogger(ProcessPool.class);
	
	private final String command;
	
	private File workDir;
	
	public ProcessPool(File workDir, String command) throws IOException {
		this.workDir = workDir;
		this.command = command;		
	}


	public Process getProcess() throws IOException, InterruptedException {
		return createProcess();
	}
	
	/**
	 * 	For now, just destroy the process. 
	 */
	public void releaseProcess(Process process) {

		if (!process.isAlive()) {
			return;
		}
		
		process.destroy();
		
		try {
			process.waitFor(WAIT_TIME_BEFORE_FORCE_DESTROY, WAIT_TIME_BEFORE_FORCE_DESTROY_UNIT);
		} catch (InterruptedException e) {
			logger.warn("interrupted while waiting process destroy");
		}
		
		if (process.isAlive()) {
			process.destroyForcibly();
		}
		
		// check again, log if still not destroyed
		try {
			process.waitFor(WAIT_TIME_BEFORE_FORCE_DESTROY, WAIT_TIME_BEFORE_FORCE_DESTROY_UNIT);
		} catch (InterruptedException e) {
			logger.warn("interrupted while waiting process force destroy");
		}
		
		if (process.isAlive()) {
			logger.warn("could not destroying process");
		}
	}
	
	private Process createProcess() throws IOException {

		ProcessBuilder builder = new ProcessBuilder(command.split(" "));
		builder.directory(workDir);
		builder.redirectErrorStream(true);
		Process p = builder.start();
		return p;
	}
	
}
