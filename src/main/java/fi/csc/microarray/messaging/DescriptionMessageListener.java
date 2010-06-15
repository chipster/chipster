package fi.csc.microarray.messaging;

import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.CountDownLatch;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.operation.OperationCategory;
import fi.csc.microarray.client.operation.OperationDefinition;
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
    
	private final CountDownLatch latch = new CountDownLatch(1);
    private Collection<OperationCategory> visibleCategories =
        new LinkedList<OperationCategory>();
    private Collection<OperationCategory> hiddenCategories =
        new LinkedList<OperationCategory>();
    private String wantedModule;
    private boolean finished = false;
    
    public DescriptionMessageListener(String wantedModule) {
        this.wantedModule = wantedModule;
    }
    
    /**
     * @return categories that are visible to end-user.
     */
    public Collection<OperationCategory> getVisibleCategories() {
        return visibleCategories;
    }

    /**
     * @return categories that are hidden from end-user, but still
     * available for execution.
     */
    public Collection<OperationCategory> getHiddenCategories() {
        return hiddenCategories;
    }
    
    public void waitForResponse() {
        try {
            latch.await();
        } catch (InterruptedException e) {
            logger.warn("interrupted while waiting for descriptions message", e);
        }
    }
    
    public void onChipsterMessage(ChipsterMessage msg) {
        if (finished) {
        	return;
        }
    	
    	ModuleDescriptionMessage descriptionMsg = (ModuleDescriptionMessage) msg;
        
    	// Only wait for one description message for the correct module
    	if (descriptionMsg.getModuleName().equals(wantedModule)) {            
        	try {
        	    parseMessage(descriptionMsg);
            } catch (ParseException e) {
                logger.warn("parsing descriptions message failed", e);
            } finally {
                finished = true;
                latch.countDown();
            }
        }
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
        
        for (Category category : categories) {
            OperationCategory op = new OperationCategory(category.getName());
            op.setColor(category.getColor());
            
            // hue, with the value .15 generates about 20 different colors
            //float hsb[] = {0f, 0.5f, 0.8f};
            //hsb[0] = (visibleCategories.size()*0.15f / 1f); 
            //op.setColor(new Color(Color.HSBtoRGB(hsb[0], hsb[1], hsb[2])));
            
            if (category.isHidden()) {
                hiddenCategories.add(op);
            } else {
                visibleCategories.add(op);
            }
            logger.debug("added category " + op.getName());
            
            // Create operation definitions for tools in this category
            for (Tool tool : category.getTools()) {
                
                SADLDescription sadl = new ChipsterSADLParser().parse(tool.getDescription());
                // logger.debug(sadl.toStringVerbose());
                
                OperationDefinition newDefinition = new OperationDefinition(sadl.getName().getID(), 
                                                                            sadl.getName().getDisplayName(), op,
                                                                            sadl.getComment(), true,
                                                                            tool.getHelpURL());
                for (Input input : sadl.inputs()) {
                    if (input.getName().isNameSet()) {
                        newDefinition.addInput(input.getName().getPrefix(), input.getName().getPostfix(), input.getType());
                    } else {
                        newDefinition.addInput(input.getName(), input.getType());
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
                // logger.debug(newDefinition.toStringVerbose());
            }
        }

        logger.debug("operation generation returning successfully");
    }
};