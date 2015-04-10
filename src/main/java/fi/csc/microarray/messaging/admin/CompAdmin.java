package fi.csc.microarray.messaging.admin;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

import javax.jms.JMSException;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.JMSMessagingEndpoint;
import fi.csc.microarray.messaging.NodeBase;
import fi.csc.microarray.messaging.admin.ServerAdminAPI.StatusReportListener;
import fi.csc.microarray.messaging.message.ServerStatusMessage;

public class CompAdmin {
	private void printHelp() {
		System.out.println("Chipster comp admin tool\n"
				+ "\n"
				+ "Usage:\n"
				+ "   java -cp path/to/jar/lib/directory/*: " + CompAdmin.class.getName() + " [-h] --config=CHIPSTER-CONFIG.XML [--wait=SECONDS] [--quiet] COMMAND\n"
				+ "\n"
				+ "Options\n"
				+ "   -h --help                     show this help and exit\n"
				+ "   --config=CHIPSTER-CONFIG.XML  Chipster manager config file\n"
				+ "   --wait=SECONDS                wait time in seconds for server responses, default 3\n"
				+ "   -q --quiet                    do not print table headers\n"
				+ "Commands\n"
				+ "   --list-comps                  list compute servers' ids, hosts and status reports\n"
				+ "   --list-jobs                   list jobs on compute servers\n"
				+ "   --cancel=JOB                  cancel a job with given id\n"		
				+ "   --stop-gracefully=COMP        wait for running jobs to complete and stop the service with given id\n"
				+ "   --stop-gracefully-slowest     gracefully stop the slowest server\n");
	}

	public int waitTime = 3; // seconds
	
	public static void main(String[] argArray) {
		
		
		LinkedList<String> args = new LinkedList<>(Arrays.asList(argArray));
		
		try {
			new CompAdmin().execute(args);
			
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}			
	}

	private CompAdminAPI compAPI;
	private JobmanagerAdminAPI jobmanagerAPI;

	private void execute(LinkedList<String> args) throws Exception {
		String config = null;
		String stopComp = null;
		String cancelJob = null;
		Integer wait = null;
		boolean listComps = false;
		boolean listJobs = false;
		boolean shutdownOne = false;
		boolean quiet = false;

		do {
			String cmd = args.poll();

			if (cmd == null) {
				cmd = "-h";
			}

			String arg = null;
			if (cmd.contains("=")) {
				String[] split = cmd.split("=");
				cmd = split[0];
				arg = split[1];
			}

			switch (cmd) {
			case "-h":
			case "--help":
				printHelp();
				System.exit(0);
				break;
			case "--config":
				config = arg;
				break;		
			case "--wait":
				wait = Integer.parseInt(arg);
				break;
			case "-q":
			case "--quiet":
				quiet = true;
				break;				
			case "--list-comps":					
				listComps = true;
				break;
			case "--list-jobs":
				listJobs = true;
				break;
			case "--cancel":
				cancelJob = arg;
				break;
			case "--stop-gracefully":
				stopComp = arg;
				break;
			case "--stop-gracefully-slowest":
				shutdownOne = true;
				break;
			default:
				break;
			}
		} while (!args.isEmpty());
		
		initConfig(config);
		
		try {

			NodeBase nodeSupport = new NodeBase() {
				public String getName() {
					return "chipster-comp-admin";
				}
			};

			ManagerConfiguration.init();
			JMSMessagingEndpoint endpoint = new JMSMessagingEndpoint(nodeSupport);					


			compAPI = new CompAdminAPI(endpoint);

			if (wait != null) {
				this.waitTime = wait;
			}

			if (listJobs) {
				listJobs(quiet);
			}
			if (listComps) {
				listComps(quiet);
			}
			if (cancelJob != null) {
				cancelJob(cancelJob);
			}
			if (stopComp != null) {
				stopGracefully(stopComp);
			}
			if (shutdownOne) {
				stopGracefully(null);
			}
			
			endpoint.close();

		} catch (MicroarrayException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (IllegalConfigurationException e) {
			e.printStackTrace();
		}
	}

	private void initConfig(String arg) throws MalformedURLException, IOException, IllegalConfigurationException {

		File config = new File(arg);
		if (!config.exists()) {			
			System.err.println("config file not found: " + config.getAbsolutePath());
			System.exit(1);
		}
		
		// server layout expects some directories to exist
		File logs = mkdir("logs");
		File security = mkdir("security");
		
		DirectoryLayout.initialiseServerLayout(Arrays.asList(new String[] {"manager"}), config.toURI().toURL().toString());
		
		if (logs != null) {
			deleteDirAndContent(logs);
		}
		
		if (security != null) {
			deleteDirAndContent(security);
		}
	}

	private void deleteDirAndContent(File dir) {
		for (File file : dir.listFiles()) {
			file.delete();
		}
		dir.delete();
	}

	private File mkdir(String name) {
		File userDir = new File(System.getProperty("user.dir"));			
		
		if (!Arrays.asList(userDir.list()).contains(name)) {
			File dir = new File(userDir.getAbsolutePath() + File.separator + name);
			if (dir.mkdir()) {
				return dir;
			}
		}
		return null;
	}

	private void stopGracefully(String compId) throws Exception {
		if (compId == null) {
			List<ServerStatusMessage> comps = new CompList().getComps(compAPI, waitTime);
			compId = comps.get(comps.size() - 1).getHostId();
		}
		
		compAPI.stopGracefullyComp(compId);
	}


	private void listComps(boolean quiet) throws Exception {
		
		if (!quiet) {
			System.out.println(ServerStatusMessage.getStringLineHeaders());
		}
		
		List<ServerStatusMessage> comps = new CompList().getComps(compAPI, waitTime);				
		
		
		for (ServerStatusMessage status : comps) {						
			System.out.println(status.toStringLine()); 						
		}		
	}
	
	class CompList {
		private List<ServerStatusMessage> comps;
		private List<ServerStatusMessage> getComps(CompAdminAPI api, int wait) throws Exception {
			comps = new ArrayList<ServerStatusMessage>();			
			api.getStatusReports(new StatusReportListener() {
				@Override
				public void statusUpdated(List<ServerStatusMessage> statuses) {					
					comps = statuses;					
				}
			}, wait);
			Thread.sleep(wait * 1000);
			return comps;
		}
	}
	
	class JobList {
		private Collection<JobsEntry> jobs;
		private Collection<JobsEntry> getJobs(JobmanagerAdminAPI jobmanagerAPI, int wait) throws Exception {
		
			jobs = jobmanagerAPI.queryRunningJobs().values();
			Thread.sleep(wait * 1000);
			return jobs;
		}
	}

	private void listJobs(boolean quiet) throws Exception {
		
		if (!quiet) {
			System.out.println(JobsEntry.getToStringHeaders());
		}
		
		Collection<JobsEntry> jobs = new JobList().getJobs(jobmanagerAPI, waitTime);
		
		for (JobsEntry job : jobs) {
			System.out.println(job.toString());						
		}	
	}
	
	private void cancelJob(String jobId) throws MicroarrayException, IOException, IllegalConfigurationException, JMSException {
		jobmanagerAPI.cancelJob(jobId);
	}
}
