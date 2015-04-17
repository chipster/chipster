package fi.csc.microarray.messaging.admin;

import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.MessagingListener;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.admin.AdminAPI.NodeStatus.Status;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.util.Strings;

/**
 * AdminAPI objects should be used only from a one thread.
 * 
 * @author Aleksi Kallio
 *
 */
public class AdminAPI {
	
	private Lock mutex = new ReentrantLock();
	
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(AdminAPI.class);

	private static final int SERVICE_CHECK_MAX_WAITTIME = 20;
	private static final int SERVICE_CHECK_MIN_WAITTIME = 5;
	
	public static final int NAME_PARAMETER_INDEX = 0;
	public static final int HOST_PARAMETER_INDEX = 1;

	public static interface AdminAPILIstener {
		public void statusUpdated(Map<String, NodeStatus> statuses);
	}
	
	public static class NodeStatus {

		public final String name;
		private HashSet<String> hosts = new HashSet<>();
		public Status status = Status.UNKNOWN;
		public int requiredCount = 0;

		public enum Status {
			UP,
			DOWN,
			UNKNOWN;
		}
		
		public NodeStatus(NodeStatus other) {
			this.name = other.name;
			this.hosts = new HashSet<>(other.hosts);
			this.status = other.status;
			this.requiredCount = other.requiredCount;
		}
		
		public NodeStatus(String name) {
			this.name = name;
		}
		
		public NodeStatus(String name, int requiredCount) {
			this(name);
			this.requiredCount = requiredCount;
		}

		public void addHost(String host) {
			hosts.add(host);
		}

		public int getCount() {
			return hosts.size();
		}

		public HashSet<String> getHosts() {
			return hosts;
		}
		
		public String toString() {
			return getClass().getSimpleName() +  " " + name + ": " + status + "(" + getCount() + " hosts)";
		}
	}

	private MessagingListener adminListener = new MessagingListener() {
		
		
		
		public void onChipsterMessage(ChipsterMessage msg) {
			mutex.lock(); 
			try {
				CommandMessage cmdMsg = (CommandMessage) msg;
				if (cmdMsg.getCommand().equals("ping-reply")) {
					String name = cmdMsg.getParameters().get(NAME_PARAMETER_INDEX);
					String host = cmdMsg.getParameters().get(HOST_PARAMETER_INDEX);
					logger.debug(name + " is up");
					NodeStatus status = nodeStatuses.get(name);
					if (status != null) {
						status.addHost(host);
						
						if (status.getCount() >= status.requiredCount) {
							status.status = Status.UP;
						}
						notifyListener();
					}		
				}
			} finally {
				mutex.unlock();
			}
		}
	};

	private MessagingTopic adminTopic;
	private Map<String, NodeStatus> nodeStatuses = new LinkedHashMap<String, NodeStatus>();
	private AdminAPILIstener listener; 
	private String errorStatus = "";
	
	public AdminAPI(MessagingTopic adminTopic, AdminAPILIstener listener) throws JMSException {
		this.listener = listener;
		this.adminTopic = adminTopic;
		adminTopic.setListener(adminListener);
		
		mutex.lock();
		try {
			this.nodeStatuses.put("authenticator", new NodeStatus("authenticator"));
			this.nodeStatuses.put("analyser", new NodeStatus("analyser"));
			this.nodeStatuses.put("filebroker", new NodeStatus("filebroker"));
			this.nodeStatuses.put("manager", new NodeStatus("manager"));
			this.nodeStatuses.put("jobmanager", new NodeStatus("jobmanager"));
			this.nodeStatuses.put("client", new NodeStatus("client"));
		} finally {
			mutex.unlock();
		}
	}

	public void setRequiredCountFor(String nodeName, int requiredCount) {
		mutex.lock();
		try {
			this.nodeStatuses.get(nodeName).requiredCount = requiredCount;
		} finally {
			mutex.unlock();
		}
	}
	
	public String getErrorStatus() {
		return errorStatus;
	}

	/**
	 * Checks that authentication service, compute service and file broker is up. 
	 * 
	 * @return true iff services are up
	 * @throws JMSException 
	 */
	public boolean areAllServicesUp(boolean fastCheck) throws InterruptedException, JMSException {
		// send ping
		ping();		
		
		// start to wait for status changes
		boolean allServicesUp = false;		
		for (int seconds = 0; ; seconds++) {
			// wait a second
			Thread.sleep(1024);
			
			// check status
			allServicesUp = areAlreadyUp();
			
			// if we are satisfied, stop waiting
			if (allServicesUp && (fastCheck || seconds >= SERVICE_CHECK_MIN_WAITTIME)) {
				break;
			}

			// if we are bored, stop waiting
			if (seconds >= SERVICE_CHECK_MAX_WAITTIME) {
				break;
			}

		}

		// update unknown statuses
		mutex.lock();
		try {
			for (NodeStatus nodeStatus : nodeStatuses.values()) {
				if (nodeStatus.status == NodeStatus.Status.UNKNOWN) {
					nodeStatus.status = NodeStatus.Status.DOWN;
				}
			}
		} finally {
			mutex.unlock();
		}
			
		return allServicesUp;
	}
	
	public String generateStatusReport() {
		String report = "";
		mutex.lock();
		try {
			for (NodeStatus status : nodeStatuses.values()) {
				report += status.name + ": count " + status.getCount() + ", host(s) " + Strings.delimit(status.getHosts(), ", ") + "\n";
			}
		} finally {
			mutex.unlock();
		}
		return report;
	}
	
	private boolean areAlreadyUp() {
		errorStatus = "";
		boolean areUp = true;

		mutex.lock();
		try {
			if (nodeStatuses.get("filebroker").status != NodeStatus.Status.UP) {
				errorStatus += " filebroker(s) not up ";
				areUp = false;
			} 

			if (nodeStatuses.get("analyser").status != NodeStatus.Status.UP) {
				errorStatus += " analyser(s) not up ";
				areUp = false;
			} 

			if (nodeStatuses.get("authenticator").status != NodeStatus.Status.UP) {
				errorStatus += " authenticator(s) not up ";
				areUp = false;
			}

			if (nodeStatuses.get("jobmanager").status != NodeStatus.Status.UP) {
				errorStatus += " jobmanager not up ";
				areUp = false;
			}
		} finally {
			mutex.unlock();
		}
		
		return areUp;
	}
	
	private void notifyListener() {
		if (listener != null) {
			mutex.lock();
			try {
				// create the copy of the map to avoid concurrent access 
				Map<String, NodeStatus> copy = new LinkedHashMap<String, NodeStatus>();
				for (String key : nodeStatuses.keySet()) {
					copy.put(key, new NodeStatus(nodeStatuses.get(key)));
				}
				listener.statusUpdated(copy);
			} finally {
				mutex.unlock();
			}
		}
	}
	
	private void ping() throws JMSException {
		// broadcast ping into admin topic
		CommandMessage pingMsg = new CommandMessage("ping");
		adminTopic.sendMessage(pingMsg);
	}	
}
