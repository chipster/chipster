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
	public final static String COMMAND_OFFER = "offer";
	public final static String COMMAND_ACCEPT_OFFER = "choose";
	public final static String COMMAND_DESCRIBE = "describe";
	public final static String COMMAND_GET_SOURCE = "get-source";

	public final static String COMMAND_URL_REQUEST ="url-request";
	public final static String COMMAND_PUBLIC_URL_REQUEST ="public-url-request";
	public final static String COMMAND_PUBLIC_FILES_REQUEST ="public-url-list-request";
	public final static String COMMAND_DISK_SPACE_REQUEST ="disk-space-request";

	
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
