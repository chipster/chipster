package fi.csc.microarray.comp;

import fi.csc.chipster.toolbox.ToolboxTool;
import fi.csc.microarray.messaging.message.GenericJobMessage;

public interface JobFactory {

	/**
	 * Creates job object using the tool description and job message, i.e., instantiates the comp
	 * job described by the description and parameterised by the parameters contained in job message. 
	 */
	public CompJob createCompJob(GenericJobMessage message, ToolboxTool tool, ResultCallback resultHandler) throws CompException;

	/**
	 * Returns true if the handler is unable to create jobs. Handler is still able to 
	 * create descriptions.
	 * 
	 * @return
	 */
	public boolean isDisabled();
}
