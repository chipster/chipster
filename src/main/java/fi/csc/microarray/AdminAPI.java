package fi.csc.microarray;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import fi.csc.microarray.AdminAPI.NodeStatus.Status;
import fi.csc.microarray.messaging.MessagingListener;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.NamiMessage;

/**
 * AdminAPI objects should be used only from a one thread.
 * 
 * @author Aleksi Kallio
 *
 */
public class AdminAPI {
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
		public String host = null;
		public Status status = Status.UNKNOWN;
		public int count = 0;
		public int requiredCount = 0;

		public enum Status {
			UP,
			DOWN,
			UNKNOWN;
		}
		
		public NodeStatus(String name) {
			this.name = name;
		}
		
		public NodeStatus(String name, int requiredCount) {
			this(name);
			this.requiredCount = requiredCount;
		}
	}

	private MessagingListener adminListener = new MessagingListener() {
		
		private Lock mutex = new ReentrantLock();
		
		public void onNamiMessage(NamiMessage msg) {
			mutex.lock(); 
			try {
				CommandMessage cmdMsg = (CommandMessage) msg;
				if (cmdMsg.getCommand().equals("ping-reply")) {
					String name = cmdMsg.getParameters().get(NAME_PARAMETER_INDEX);
					String host = cmdMsg.getParameters().get(HOST_PARAMETER_INDEX);
					logger.debug(name + " is up");
					NodeStatus status = nodeStatuses.get(name);
					if (status != null) {
						if (status.host == null) {
							status.host = host;
						} else {
							status.host += ", " + host;
						}
						status.count += 1;
						if (status.count >= status.requiredCount) {
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
		this.nodeStatuses.put("authenticator", new NodeStatus("authenticator"));
		this.nodeStatuses.put("analyser", new NodeStatus("analyser"));
		this.nodeStatuses.put("client", new NodeStatus("client"));
	}

	public void setRequiredCountFor(String nodeName, int requiredCount) {
		this.nodeStatuses.get(nodeName).requiredCount = requiredCount;
	}
	
	public String getErrorStatus() {
		return errorStatus;
	}

	/**
	 * Checks that analyser and frontend are up.
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
		for (NodeStatus nodeStatus : nodeStatuses.values()) {
			if (nodeStatus.status == NodeStatus.Status.UNKNOWN) {
				nodeStatus.status = NodeStatus.Status.DOWN;
			}
		}
			
		return allServicesUp;
	}
	
	public String generateStatusReport() {
		String report = "";
		for (NodeStatus status : nodeStatuses.values()) {
			report += status.name + ": count " + status.count + ", host(s) " + status.host + "\n";
		}
		return report;
	}
	
	private boolean areAlreadyUp() {
		errorStatus = "";
		boolean areUp = true;
		
		if (nodeStatuses.get("analyser").status != NodeStatus.Status.UP) {
			errorStatus += " analyser(s) not up ";
			areUp = false;
		} 
		
		if (nodeStatuses.get("authenticator").status != NodeStatus.Status.UP) {
			errorStatus += " authenticator(s) not up ";
			areUp = false;
		}
		
		return areUp;
	}
	
	private void notifyListener() {
		if (listener != null) {
			listener.statusUpdated(nodeStatuses);
		}
	}
	
	private void ping() throws JMSException {
		// broadcast ping into admin topic
		CommandMessage pingMsg = new CommandMessage("ping");
		adminTopic.sendMessage(pingMsg);
	}	
}
