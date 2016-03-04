package fi.csc.microarray.messaging.message;

import java.util.List;
import java.util.Set;
import java.util.UUID;

import fi.csc.microarray.comp.ToolDescription;
import fi.csc.microarray.messaging.message.JobMessage.ParameterSecurityPolicy;
import fi.csc.microarray.messaging.message.JobMessage.ParameterValidityException;

/**
 * Use this interface instead of the original JobMessage on the comp. The
 * original class has a JMS dependency, which we don't want in the new system
 * which is communicating with REST APIs.
 * 
 * @author klemela
 */
public interface GenericJobMessage {

	public String getJobId();

	public String getUsername();

	public String getToolId();

	public Set<String> getKeys() throws Exception;

	public String getId(String fileName);

	public String getName(String fileName);

	public List<String> getParameters(ParameterSecurityPolicy securityPolicy, ToolDescription description) throws ParameterValidityException;

	public UUID getSessionId();

}
