package fi.csc.microarray.comp;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.Date;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.LinkedBlockingQueue;
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

	static final Logger logger = Logger.getLogger(ProcessPool.class);
	
	private BlockingQueue<NamiProcess> availableProcesses;
	private ConcurrentMap<Integer, NamiProcess> inUseProcesses;
	
	private final int poolSizeMin;
	private final int poolSizeMax;
	private final int poolTimeout;
	private final int processUseCountMax;
	private final int processLifetimeMax;
	private final String rCommand;
	
	private File workDir;
	
	private static final String RECYCLE_SUCCESFUL_STRING = "recycling-succesful";
	
	
	/**
	 * 
	 * Wrapper class for storing metadata about a Process.
	 * 
	 */
	private class NamiProcess {
		
		private Process process;
		private int useCount;
		private Date creationTime;
		
		public NamiProcess(Process p) {
			this.process = p;
			useCount = 0;
			creationTime = new Date(System.currentTimeMillis());
		}
		
		public Process getProcess() {
			return this.process;
		}

		public int getUseCount() {
			return this.useCount;
		}
		
		public Date getCreationTime() {
			return this.creationTime;
		}
		
		public void increaseUseCount() {
			this.useCount++;
		}
		
	}
	
	
	
	public ProcessPool(File workDir, String command, int poolSizeMin, int poolSizeMax, int poolTimeout, int processUseCountMax, int processLifetimeMax) throws IOException {
		this.workDir = workDir;
		this.rCommand = command;
		this.poolSizeMin = poolSizeMin;
		this.poolSizeMax = poolSizeMax;
		this.poolTimeout = poolTimeout;
		this.processUseCountMax = processUseCountMax;
		this.processLifetimeMax = processLifetimeMax;

		// initialize pool structures
		this.availableProcesses = new LinkedBlockingQueue<NamiProcess>(); 
		this.inUseProcesses = new ConcurrentHashMap<Integer, NamiProcess>();
	
		// populate pool with new processes  
		for (int i = 0; i < poolSizeMin; i++) {
			availableProcesses.add(createProcess());
		}
		
		logger.debug("R process pool initialized, processes available: " + availableProcesses.size());
		
	}


	public Process getProcess() throws IOException, InterruptedException {
		NamiProcess nProcess;
		
		// try to get a process
		nProcess = availableProcesses.poll();
		
		// no processes were available, try to create a new one
		if (nProcess == null) {
			
			// still room for more processes, create one
			if (availableProcesses.size() + inUseProcesses.size() < poolSizeMax) {
				nProcess = createProcess();
			}

			// already max number of processes created, wait for one to be available
			else {
				nProcess = availableProcesses.poll(poolTimeout, TimeUnit.SECONDS);
				
				// finally got the process or timeout?
				if (nProcess == null) {
					throw new IOException("Timeout when getting an R process.");
				}
			}
		}
		
		// the process is now in use
		nProcess.increaseUseCount();
		inUseProcesses.put(nProcess.getProcess().hashCode(), nProcess);
		return nProcess.getProcess();
	}
	
	/**
	 * 
	 * @param process
	 * @param recycle true if the process should be recycled, false is used for
	 * processes which are known to be dead or having problems
	 * @throws IOException
	 */
	public void releaseProcess(Process process, boolean recycle) throws IOException {

		
		// make sure the process originated in this pool
		NamiProcess nProcess = inUseProcesses.get(process.hashCode());
		if (nProcess == null) {
			throw new IOException("Trying to release an unknown process.");
		}
		
		
		// check the process
		boolean processOk = true;
		boolean processAlive = false;
		
		// no recycling 
		if (!recycle) {
			logger.debug("Process " + nProcess.getProcess().hashCode() + " not recycled as requested.");
		} 
		// check process use count
		else if (nProcess.getUseCount() >= processUseCountMax) {
			processOk = false;
			logger.debug("Process " + nProcess.getProcess().hashCode() + " has been used for " + nProcess.getUseCount() + " times and is therefore not recycled.");
		} 
		// check process lifetime
		else if ((System.currentTimeMillis() - nProcess.getCreationTime().getTime()) >= processLifetimeMax * 1000) {
			processOk = false;
			logger.debug("Process " + nProcess.getProcess().hashCode() + " has been running over " + processLifetimeMax + " seconds and is therefore not recycled.");
		}
		// check if process is alive
		else {
			try {
				process.exitValue();
			} catch (IllegalThreadStateException itse) {
				// if we got here, process is still running
				processAlive = true;
			}
			if (!processAlive) {
				logger.debug("Process " + nProcess.getProcess().hashCode() + " is dead and is therefore not recycled.");
			}
		}

		
		// recycle the process or replace it with a new one 
		if (recycle && processOk && processAlive) {
			logger.debug("Recycling process " + nProcess.getProcess().hashCode() + ".");

			// make sure that STDOUT and STDERR are clean

			
			CountDownLatch recycleLatch = new CountDownLatch(1);
			ProcessMonitor recycleMonitor = new ProcessMonitor(process, recycleLatch);
			new Thread(recycleMonitor).start();
			
			BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(process.getOutputStream()));
			// TODO add garbage collecting here?
			writer.write("rm(list=objects())");
			writer.newLine();
			writer.write("setwd(\"" + workDir.getAbsolutePath() + "\")");
			writer.newLine();
			writer.write("print(\"" + RECYCLE_SUCCESFUL_STRING + "\")");
			writer.newLine();
			writer.flush();
			
			try {
				recycleLatch.await(10, TimeUnit.SECONDS);
			} catch (InterruptedException e) {
				processOk = false;
			}

			if (!processOk || recycleLatch.getCount() > 0 || !recycleMonitor.processOk() || process.getErrorStream().available() > 0) {
				process.destroy();
				nProcess = createProcess();
			}
			
		
		} else {
			process.destroy();
			nProcess = createProcess();
		}

		
		// make the (possibly new) process available again, if there is room for it,
		if (availableProcesses.size() < poolSizeMin && availableProcesses.size() + inUseProcesses.size() <= poolSizeMax) {
			availableProcesses.add(nProcess);
		} 
		// otherwise destroy it
		else {
			nProcess.getProcess().destroy();
		}
		
		// remove the possibly recycled process from inUse
		inUseProcesses.remove(process.hashCode());
		
		logger.debug("Available processes: " + availableProcesses.size() + ", in use: " + inUseProcesses.size());

	
	}
	

	private NamiProcess createProcess() throws IOException {

		logger.debug("Creating a new R process.");
		ProcessBuilder builder = new ProcessBuilder(rCommand.split(" "));
		//Process p = Runtime.getRuntime().exec(rCommand, null, workDir);
		builder.directory(workDir);
		builder.redirectErrorStream(true);
		Process p = builder.start();
		return new NamiProcess(p);
	}
	
	
	private class ProcessMonitor implements Runnable {

		private CountDownLatch latch;
		private Process process;
		private boolean processOk = false;
		
		public ProcessMonitor(Process process, CountDownLatch latch) {
			this.process = process;
			this.latch = latch;
		}
		
		public boolean processOk() {
			return processOk;
		}
		
		
		public void run() {
			BufferedReader reader = new BufferedReader(new InputStreamReader(this.process.getInputStream()));
			
			boolean readMore = true;
			try {
				for (String line = reader.readLine(); readMore ; line = reader.readLine()) {
					
					// read end of stream --> error
					if (line == null ) {
						processOk = false;
						readMore = false;
					} 
					
					// read recycle successful
					else if (line.contains(RECYCLE_SUCCESFUL_STRING)) {
						processOk = true;
						readMore = false;
					}
					
					// read normal output
					else {
						
					}
				}
			} catch (IOException e) {
				processOk = false;
			}

			this.latch.countDown();
		}
	
	}

	
	
	
	
	
	
	
	
}
