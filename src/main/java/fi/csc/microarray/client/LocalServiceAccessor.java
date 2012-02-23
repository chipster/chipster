package fi.csc.microarray.client;

import java.util.Collection;
import java.util.LinkedList;

import javax.jms.JMSException;

import fi.csc.chipster.tools.ngs.LocalNGSPreprocess;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.ToolCategory;
import fi.csc.microarray.client.operation.ToolModule;
import fi.csc.microarray.client.tasks.LocalTaskExecutor;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.filebroker.SimpleFileBrokerClient;
import fi.csc.microarray.messaging.SourceMessageListener;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.messaging.message.FeedbackMessage;
import fi.csc.microarray.module.Module;
import fi.csc.microarray.module.chipster.ChipsterSADLParser;

public class LocalServiceAccessor implements ServiceAccessor {

	private static final String LOCAL_CATEGORY_NAME = "Local tools";
	private static final String LOCAL_MODULE_NAME = "local";

	private DataManager manager;

	private LinkedList<ToolModule> modules = new LinkedList<ToolModule>();

	
	
	@Override
	public String checkRemoteServices() throws Exception {
		return ServiceAccessor.ALL_SERVICES_OK;
	}

	@Override
	public void close() throws Exception {
		// do nothing
	}


	/**
	 * FIXME get tool list from configs
	 * FIXME put the code for creating the OperationDefinition to a place where both this and
	 * DescriptionMessageListener can use it
	 * 
	 */
	@Override
	public void fetchDescriptions(Module primaryModule) throws Exception {
		ToolModule module = new ToolModule(LOCAL_MODULE_NAME);
		this.modules.add(module);
		
		// for each tool, parse the SADL and create the OperationDefinition
        SADLDescription sadl = new ChipsterSADLParser().parse(LocalNGSPreprocess.getSADL());
		ToolCategory category = new ToolCategory(LOCAL_CATEGORY_NAME);
		module.addHiddenToolCategory(category);
		
        OperationDefinition od = new OperationDefinition(sadl.getName().getID(), 
                                                                    sadl.getName().getDisplayName(), category,
                                                                    sadl.getComment(), true,
                                                                    null);
        for (Input input : sadl.inputs()) {
            if (input.getName().isNameSet()) {
                od.addInput(input.getName().getPrefix(), input.getName().getPostfix(), input.getName().getDisplayName(), input.getComment(), input.getType(), input.isOptional());
            } else {
                od.addInput(input.getName(), input.getComment(), input.getType(), input.isOptional());
            }
        }

        od.setOutputCount(sadl.outputs().size());
        for (Parameter parameter : sadl.parameters()) {
            od.addParameter(fi.csc.microarray.client.
                                       operation.parameter.Parameter.createInstance(
                parameter.getName(), parameter.getType(), parameter.getSelectionOptions(),
                parameter.getComment(), parameter.getFrom(), parameter.getTo(),
                parameter.getDefaultValues(), parameter.isOptional()));      
        }
	}

	
	@Override
	public FileBrokerClient getFileBrokerClient() throws Exception {
		return new SimpleFileBrokerClient();
	}

	
	@Override
	public TaskExecutor getTaskExecutor() {
		try {
			return new LocalTaskExecutor(manager);
		} catch (JMSException e) {
			throw new RuntimeException(e);
		}
	}

	@Override
	public void initialise(DataManager manager, AuthenticationRequestListener authenticationRequestListener) throws Exception {
		this.manager = manager;
		// we are not interested in authenticationRequestListener
	}

	@Override
	public SourceMessageListener retrieveSourceCode(String id) throws Exception {
		throw new UnsupportedOperationException("not supported in standalone mode");
	}

	@Override
	public void sendFeedbackMessage(FeedbackMessage message) throws Exception {
		throw new UnsupportedOperationException("not supported in standalone mode");
	}

	@Override
	public Collection<ToolModule> getModules() {
		return modules;
	}

	@Override
	public boolean isStandalone() {
		return true;
	}

}
