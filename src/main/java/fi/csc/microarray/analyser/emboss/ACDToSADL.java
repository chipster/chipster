package fi.csc.microarray.analyser.emboss;

import java.util.LinkedList;

import fi.csc.microarray.description.ParsedVVSADL;

/**
 * Conversion of ACD abstraction to SADL abstraction.
 * 
 * @author naktinis
 *
 */
public class ACDToSADL {
	
	private ACDDescription acd;
	
	public ACDToSADL(ACDDescription acd) {
		this.acd = acd;
	}

	/**
	 * Analyse a given ACD object and store it as a SADL abstraction.
	 * 
	 * @return SADL object.
	 */
	public ParsedVVSADL convert() {
	    ParsedVVSADL sadl = new ParsedVVSADL(acd.getName(), acd.getGroups().get(0),
	                                         acd.getDescription());
	    
	    // Get all input parameters
	    // We are also safe from trying to include non-input parameters in input
	    //     section (such as toggle), since SADLParameterCreator does type-checking.
	    // TODO: deal with non-input parameters in input and output sections
	    
	    LinkedList<ACDParameter> params = acd.getParameters();
	    for (ACDParameter param : params) {
	        SADLParameterCreator.createAndAdd(param, sadl);
        }
	    
	    
	    /*
	    LinkedList<ACDParameter> inputParams = acd.getParameters("input", null, true);
	    for (ACDParameter inputParam : inputParams) {
	        ParsedVVSADL.Input sadlInput = SADLParameterCreator.createInput(inputParam);
            sadl.addInput(sadlInput);
        }

	    // Get all simple parameters and store them
	    // NOTE: currently does not fetch advanced parameters
        LinkedList<ACDParameter> displayedParams = acd.getParameters("required", null, true);
        displayedParams.addAll(acd.getParameters("additional", null, true));
        for (ACDParameter param : displayedParams) {
            ParsedVVSADL.Parameter sadlParam = SADLParameterCreator.createParameter(param);
            sadl.addParameter(sadlParam);
        }

        // Get all output parameters
        LinkedList<ACDParameter> outputParams = acd.getParameters("output", null, true);
        for (ACDParameter outputParam : outputParams) {
            String sadlOutput = SADLParameterCreator.createOutput(outputParam);
            sadl.addOutput(sadlOutput);
        }*/
        
	    return sadl;
    }
}