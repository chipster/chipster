package fi.csc.microarray.comp;

import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Timer;
import java.util.TimerTask;

import org.apache.log4j.Logger;

import fi.csc.microarray.comp.ProcessUtils.ProcessResourceMonitor;

/**
 * <h1>Monitor resource usage of unix-like processes</h1>
 * 
 * <p>Get the list of current processes from the ProcessProvider and collect the 
 * combined memory usage of each process and its child processes every second. When the process disappears
 * from the list of ProcessProvider, it's monitoring results are removed as well.</p>
 * 
 * <p>Java doesn't provide APIs for doing this, so external tools are used instead. There are multiple things
 * that can go wrong or cause inaccuracies, so be prepared for null results.</p>
 * <ul>
 * <li> pid of the process is figured out with the Java Reflection API</li>
 * <li> child processes are parsed from the output of the pgrep command</li>
 * <li> memory usage is parsed from the output of the ps command</li>
 * <li> tracking processes and memory usage by sampling is inaccurate</li>
 * </ul>
 * 
 * <p>Tested only on Linux and OSX.</p>
 * 
 * @author klemela
 *
 */
public class ResourceMonitor {
	
	public static interface ProcessProvider {
		public Collection<Process> getRunningJobProcesses();
	}

	static final Logger logger = Logger.getLogger(ResourceMonitor.class);
	
	private HashMap<Process, ProcessResourceMonitor> monitors = new HashMap<>();	
	private Timer resourceMonitorTimer;
	
	private ProcessProvider processProvider;

	public ResourceMonitor(ProcessProvider processProvider, int monitoringInterval) {
		if (monitoringInterval >= 0) {
			this.processProvider = processProvider;
			
			resourceMonitorTimer = new Timer(true);
			resourceMonitorTimer.schedule(new ResourceMonitorTask(), monitoringInterval, monitoringInterval);
		}
	}
	
	public class ResourceMonitorTask extends TimerTask {
		

		@Override
		public void run() {
			try {
				
				long t = System.currentTimeMillis();
				
				Collection<Process> runningProcesses = processProvider.getRunningJobProcesses();
				
				logger.debug("running processes " + runningProcesses.size());
				
				// remove monitor if the job isn't running anymore
				monitors.keySet().retainAll(runningProcesses);
				
				// create or update monitor of each running job
				for (Process process : runningProcesses) {
					if (!monitors.containsKey(process)) {
						monitors.put(process, new ProcessResourceMonitor(process));					
					}						
					monitors.get(process).update();
				}
				
				long dt = (System.currentTimeMillis() - t);
				if (dt > 500) {
					// consider getting information of all pids with a single ps process if this happens often
					logger.warn("process monitoring took " + (System.currentTimeMillis() - t) + "ms");
				}
				
			} catch (IOException e) {
				logger.error("failed to monitor process resource usage", e);
			}
		}	
	}
	
	public Long getMaxMem(Process process) {
		// return null if monitoring is disabled (this.monitors is still initialized) 
		ProcessResourceMonitor monitor = monitors.get(process);
		if (monitor == null) {
			return null;
		}
		return monitor.getMaxMem();
	}
}
