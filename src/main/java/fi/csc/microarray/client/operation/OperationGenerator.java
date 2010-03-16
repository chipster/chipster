package fi.csc.microarray.client.operation;

import java.awt.Color;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLDescription.Parameter;

public class OperationGenerator {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(OperationGenerator.class);

	public Map<String, OperationCategory> generate(List<SADLDescription> descriptions) {
		logger.debug("generating operations from " + descriptions.size() + " descriptions");
		LinkedHashMap<String, OperationCategory> parsedCategories = new LinkedHashMap<String, OperationCategory>();
		for (SADLDescription description : descriptions) {
			logger.debug("processing " + description.getName());
			
			// if the category doesn't exist yet, create it
			if (parsedCategories.get(description.getCategory()) == null) {
				OperationCategory op = new OperationCategory(description.getCategory());
				
				// predefined colors are used normally ...
				if (parsedCategories.size() + 1  < VisualConstants.CATEGORY_COLORS.length) {
					op.setColor(VisualConstants.CATEGORY_COLORS[parsedCategories.size() + 1]);
					
				
				} else {
					// ... and generated colors when there are not enough predefined ones
					
					// decide the color of category for the visualization purposes			
					float hsb[] = {0f, 0.5f, 0.8f};
					
					// hue, with the value .15 generates about 20 different colors
					hsb[0] = (parsedCategories.size()*0.15f / 1f); 
					
					op.setColor(new Color(Color.HSBtoRGB(hsb[0], hsb[1], hsb[2])));
				}
				
				parsedCategories.put(op.getName(), op);
				logger.debug("added category " + op.getName());
			}
			
			// create the actual definition
			OperationDefinition newDefinition = new OperationDefinition(description.getName().getID(), parsedCategories.get(description.getCategory()), description.getComment(), true);
			logger.debug("added operation " + newDefinition.getName() + " to " + newDefinition.getCategoryName());
			
			for (Input input : description.inputs()) {
				if (input.getName().isNameSet()) {
					newDefinition.addInput(input.getName().getPrefix(), input.getName().getPostfix(), input.getType());
				} else {
					newDefinition.addInput(input.getName().getID(), input.getType());
				}
			}
			
			logger.debug("added " + description.inputs().size() + " inputs");

			newDefinition.setOutputCount(description.outputs().size());
			for (Parameter parameter : description.parameters()) {
				newDefinition.addParameter(fi.csc.microarray.client.
				                           operation.parameter.Parameter.createInstance(
			        parameter.getName().getID(), parameter.getType(), parameter.getSelectionOptions(),
			        parameter.getComment(), parameter.getFrom(), parameter.getTo(),
			        parameter.getDefaultValue()));		
			}
		}
		
		logger.debug("operation generation returning successfully");
		return parsedCategories;
	}
}
