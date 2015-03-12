package fi.csc.microarray.messaging.message;

import java.util.List;

import javax.jms.JMSException;
import javax.jms.MapMessage;

/**
 * For sending a named command and a list of parameters.
 *  
 * @author akallio
 */
public class CommandMessage extends ParameterMessage {

	
	private final static String KEY_COMMAND = "command";

	public final static String COMMAND_CANCEL = "cancel";
	public final static String COMMAND_ACK = "acknowledge";
	public final static String COMMAND_EOT = "end-of-transmission";
	public final static String COMMAND_OFFER = "offer";
	public final static String COMMAND_ACCEPT_OFFER = "choose";
	public final static String COMMAND_DESCRIBE = "describe";
	public final static String COMMAND_GET_SOURCE = "get-source";

	public final static String COMMAND_NEW_URL_REQUEST ="new-url-request";
	public final static String COMMAND_GET_URL ="get-url";

	@Deprecated
	public final static String COMMAND_PUBLIC_URL_REQUEST ="public-url-request";
	public final static String COMMAND_PUBLIC_FILES_REQUEST ="public-url-list-request";
	public final static String COMMAND_DISK_SPACE_REQUEST ="disk-space-request";

	public static final String COMMAND_MOVE_FROM_CACHE_TO_STORAGE = "move-from-cache-to-storage";
	public static final String COMMAND_IS_AVAILABLE = "is-available";
	public static final String COMMAND_STORE_SESSION = "store-session";
	public static final String COMMAND_REMOVE_SESSION = "remove-session";
	public static final String COMMAND_LIST_SESSIONS = "list-sessions";
	public static final String COMMAND_GET_STATUS_REPORT = "get-status-report";

	public static final String COMMAND_FILE_OPERATION_SUCCESSFUL = "file-operation-successful";
	public static final String COMMAND_FILE_OPERATION_FAILED = "file-operation-failed";
	public static final String COMMAND_FILE_OPERATION_DENIED = "file-operation-denied";

	public static final String COMMAND_LIST_STORAGE_USAGE_OF_USERS = "get-storage-usage-by-user";
	public static final String COMMAND_LIST_STORAGE_USAGE_OF_SESSIONS = "get-sessions-for-user";
	public static final String COMMAND_GET_STORAGE_USAGE_TOTALS = "get-sessions-for-session-name";
	
	public static final String COMMAND_LIST_RUNNING_JOBS = "get-running-jobs";
	public static final String COMMAND_GET_COMP_STATUS = "get-comp-status";
	public static final String COMMAND_STOP_GRACEFULLY_COMP = "shutdown-comp";

	public static final String COMMAND_GET_JOB = "get-job";

	public static final String COMMAND_PURGE_OLD_JOBS = "purge-old-jobs";	
	
	private String command;
	

	/**
	 * When only one kind of commands are passed to a given topic, this constructor
	 * can be used to create anonymous command.
	 */
	public CommandMessage() {
		this(null, null);
	}
	
	public CommandMessage(String command) {
		this(command, null);		
	}
	
	public CommandMessage(String command, List<String> parameters) {
		super(parameters);
		this.command = command;		
	}
	
	@Override
	public void unmarshal(MapMessage from) throws JMSException {
		super.unmarshal(from);
		this.command = from.getString(KEY_COMMAND);
	}
	
	@Override
	public void marshal(MapMessage to) throws JMSException {
		super.marshal(to);
		to.setString(KEY_COMMAND, command);
	}
	
	/**
	 * Returns the command or null if it is anonymous.
	 */
	public String getCommand() {
		return command;
	}
	
	/**
	 * Gets parameters in the order they were inserted.
	 * No safety precautions, should only be used in messaging
	 * between system components (never for user input).
	 */
	public List<String> getParameters() {
		return super.getParameters();
	}
}
