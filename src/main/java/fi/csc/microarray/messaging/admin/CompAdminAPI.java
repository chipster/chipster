package fi.csc.microarray.messaging.admin;

import java.io.IOException;
import java.util.Collection;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.SuccessMessageListener;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.ParameterMessage;
import fi.csc.microarray.messaging.message.SuccessMessage;

/**
 * This class uses JMS messages to send data queries and converts result messages to
 * Java objects. The methods wait for the results, turning asynchronous messages to 
 * blocking method calls.
 * 
 * @author klemela
 */
public class CompAdminAPI extends ServerAdminAPI {
		
	private static final Logger logger = Logger.getLogger(CompAdminAPI.class);
	
	public static interface JobsListener {
		public void statusUpdated(Collection<JobsEntry> collection);
	}

	public CompAdminAPI() throws IOException, IllegalConfigurationException, MicroarrayException, JMSException {
		super(Topics.Name.COMP_ADMIN_TOPIC, "comp-admin");
	}
	
	public void stopGracefullyComp(String compId) throws MicroarrayException {
		SuccessMessageListener replyListener = new SuccessMessageListener();  
				
		try {
			CommandMessage removeRequestMessage = new CommandMessage(CommandMessage.COMMAND_STOP_GRACEFULLY_COMP);
			removeRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_HOST_ID, compId); 
			getTopic().sendReplyableMessage(removeRequestMessage, replyListener);

			SuccessMessage reply = replyListener.waitForReply(TIMEOUT, TIMEOUT_UNIT);
			
			checkSuccessMessage(reply, "stop comp gracefully");
												
		} catch (JMSException e) {
			logger.error("stopping comp gracefully failed", e);
		}
	}
}