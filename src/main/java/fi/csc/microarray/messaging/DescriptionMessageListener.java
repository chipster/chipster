package fi.csc.microarray.messaging;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.concurrent.CountDownLatch;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.ToolCategory;
import fi.csc.microarray.client.operation.ToolModule;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.ModuleDescriptionMessage;
import fi.csc.microarray.messaging.message.ModuleDescriptionMessage.Category;
import fi.csc.microarray.messaging.message.ModuleDescriptionMessage.Tool;
import fi.csc.microarray.module.chipster.ChipsterSADLParser;

public class DescriptionMessageListener extends TempTopicMessagingListenerBase {
    
	private static final Logger logger = Logger.getLogger(DescriptionMessageListener.class);
    
	private HashMap<String, ToolModule> modules = new LinkedHashMap<String, ToolModule>();

	private final CountDownLatch latch = new CountDownLatch(1);
    private String[] requiredModules;
    private boolean finished = false;

    private boolean[] isModuleLoaded;
	private ModuleDescriptionMessage[] messages;
    
    public DescriptionMessageListener(String[] requiredModules) {
        this.requiredModules = requiredModules;
        this.isModuleLoaded = new boolean[requiredModules.length];
        this.messages = new ModuleDescriptionMessage[requiredModules.length];
        Arrays.fill(this.isModuleLoaded, false);
    }
    

	public Collection<ToolModule> getModules() {
		return modules.values();
	}

    public void waitForResponse() {
        try {
            latch.await();
        } catch (InterruptedException e) {
            logger.warn("interrupted while waiting for descriptions message", e);
        }
    }
    
    public void onChipsterMessage(ChipsterMessage cMsg) {
    	if (finished) {
    		return;
    	}

    	ModuleDescriptionMessage descriptionMsg = (ModuleDescriptionMessage) cMsg;

    	// Filter those modules we are interested in
    	for (int i = 0; i < requiredModules.length; i++) {
    		if (requiredModules[i].equals(descriptionMsg.getModuleName()) && isModuleLoaded[i] == false) {

    			// Store the descriptions
    			isModuleLoaded[i] = true; // mark it loaded if we have at least tried to load it
    			messages[i] = descriptionMsg;
    			
    			// Are we done now?
    			if (isAllModulesLoaded()) {

    				for (ModuleDescriptionMessage mdMsg : messages) {
    	    			try {
    	    				parseMessage(mdMsg);
    	    			} catch (ParseException e) {
    	    				logger.warn("parsing descriptions message failed", e);
    	    			}
    				}
    				finished = true;
    				latch.countDown();
    			}
    		}
    	}
    }
    
    private boolean isAllModulesLoaded() {
    	for (boolean isLoaded : isModuleLoaded) {
    		if (!isLoaded) {
    			return false;
    		}
    	}
    	return true;
    }

    /**
     * Prepare operation category lists.
     * 
     * @param descriptionMsg
     * @throws ParseException
     */
    private void parseMessage(ModuleDescriptionMessage descriptionMsg) throws ParseException {
        
        // Fetch descriptions from the ModuleDescriptionMessage
        List<Category> categories = descriptionMsg.getCategories();
        logger.debug("generating operations from " + categories.size() + " categories");
        ToolModule module = new ToolModule(descriptionMsg.getModuleName());
        
        for (Category category : categories) {
            ToolCategory toolCategory = new ToolCategory(category.getName());
            toolCategory.setColor(category.getColor());
            toolCategory.setModule(module);
            
            if (category.isHidden()) {
                module.addHiddenToolCategory(toolCategory);
            } else {
                module.addVisibleToolCategory(toolCategory);
            }
            logger.debug("added category " + toolCategory.getName());
            
            // Create operation definitions for tools in this category
            for (Tool tool : category.getTools()) {
                
                SADLDescription sadl = new ChipsterSADLParser().parse(tool.getDescription());
                
                // check for duplicate id in this module
                if (module.getOperationDefinition(sadl.getName().getID()) != null) {
                	logger.warn("ignoring tool with duplicate id: " + sadl.getName().getID());
                	Session.getSession().getApplication().showDialog("Ignoring tool with duplicate id " + sadl.getName().getID() , 
                				"This message is meant for the Chipster service administrator, who should check the tool configurations. Don't worry about it.", "", Severity.INFO, false);
                	continue;
                }
                
                
                OperationDefinition newDefinition = new OperationDefinition(sadl.getName().getID(), 
                                                                            sadl.getName().getDisplayName(), toolCategory,
                                                                            sadl.getComment(), true,
                                                                            tool.getHelpURL());
                for (Input input : sadl.inputs()) {
                    if (input.getName().isNameSet()) {
                        newDefinition.addInput(input.getName().getPrefix(), input.getName().getPostfix(), input.getName().getDisplayName(), input.getComment(), input.getType(), input.isOptional());
                    } else {
                        newDefinition.addInput(input.getName(), input.getComment(), input.getType(), input.isOptional());
                    }
                }
    
                newDefinition.setOutputCount(sadl.outputs().size());
                for (Parameter parameter : sadl.parameters()) {
                    newDefinition.addParameter(fi.csc.microarray.client.
                                               operation.parameter.Parameter.createInstance(
                        parameter.getName(), parameter.getType(), parameter.getSelectionOptions(),
                        parameter.getComment(), parameter.getFrom(), parameter.getTo(),
                        parameter.getDefaultValues(), parameter.isOptional()));      
                }
            }
        }
        
        modules.put(module.getModuleName(), module);

        logger.debug("operation generation returning successfully");
    }
};