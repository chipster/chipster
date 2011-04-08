package fi.csc.microarray.client.operation;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.messaging.message.ModuleDescriptionMessage;
import fi.csc.microarray.messaging.message.ModuleDescriptionMessage.Category;
import fi.csc.microarray.messaging.message.ModuleDescriptionMessage.Tool;
import fi.csc.microarray.module.chipster.ChipsterSADLParser;

public class OperationGenerator {
    
	// Logger for this class
	private static final Logger logger = Logger.getLogger(OperationGenerator.class);
	
	public Map<String, OperationCategory> generateFromMessage(ModuleDescriptionMessage descriptionMsg) throws ParseException {
        LinkedHashMap<String, OperationCategory> parsedCategories =
                            new LinkedHashMap<String, OperationCategory>();
	    // Fetch descriptions from the DescriptionMessage
	    List<Category> categories = descriptionMsg.getCategories();
	    
	    logger.debug("generating operations from " + categories.size() + " descriptions");
	    
	    for (Category category : categories) {
	        OperationCategory op = new OperationCategory(category.getName());
	        op.setColor(category.getColor());
            
            // hue, with the value .15 generates about 20 different colors
            //float hsb[] = {0f, 0.5f, 0.8f};
            //hsb[0] = (parsedCategories.size()*0.15f / 1f); 
            //op.setColor(new Color(Color.HSBtoRGB(hsb[0], hsb[1], hsb[2])));
	        
            parsedCategories.put(op.getName(), op);
            logger.debug("added category " + op.getName());
            
            // Create operation definitions for tools in this category
            for (Tool tool : category.getTools()) {
                
                SADLDescription sadl = new ChipsterSADLParser().parse(tool.getDescription());
                
                OperationDefinition newDefinition = new OperationDefinition(sadl.getName().getID(), 
                															sadl.getName().getDisplayName(), op,
                                                                            sadl.getComment(), true,
                                                                            tool.getHelpURL());
                for (Input input : sadl.inputs()) {
                    if (input.getName().isNameSet()) {
                        newDefinition.addInput(input.getName().getPrefix(), input.getName().getPostfix(), input.getName().getDisplayName(), input.getComment(), input.getType());
                    } else {
                        newDefinition.addInput(input.getName(), input.getComment(), input.getType());
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
                logger.debug(newDefinition.toStringVerbose());
            }
	    }

	    logger.debug("operation generation returning successfully");
	    return parsedCategories;
	}
}
