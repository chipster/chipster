package fi.csc.microarray.analyser.emboss;

import java.util.LinkedList;

import fi.csc.microarray.description.SADLDescription;

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
	public SADLDescription convert() {
        SADLDescription sadl = new SADLDescription(acd.getName(), acd.getGroups().get(0),
	                                               acd.getDescription());
	    
	    // Get all input parameters
	    // We are also safe from trying to include non-input parameters in input
	    //     section (such as toggle), since SADLParameterCreator does type-checking.
	    // TODO: deal with non-input parameters in input and output sections
	    
	    LinkedList<ACDParameter> params = acd.getParameters();
	    for (ACDParameter param : params) {
	        SADLParameterCreator.createAndAdd(param, sadl);
        }
        
	    return sadl;
    }
}