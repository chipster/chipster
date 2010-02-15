package fi.csc.microarray.analyser.emboss;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedHashMap;

import org.apache.regexp.RE;
import org.emboss.jemboss.parser.AcdFunResolve;
import org.emboss.jemboss.parser.ParseAcd;

import fi.csc.microarray.description.ParsedVVSADL;

/**
 * Manipulation of ACD files which describe EMBOSS applications.
 * 
 * @author naktinis
 *
 */
public class ACDReader {
	
	private String filename;
	private HashMap <String, String>appAttrs;
	
	public ACDReader (String filename) {
		this.filename = filename;
		
		// Application-level attributes that we want to extract
		appAttrs = new HashMap<String, String>();
		appAttrs.put("application", "");
		appAttrs.put("documentation", "");
		appAttrs.put("groups", "");
	}

	/**
	 * Analyse a given ACD file and store it in some internal
	 * format.
	 *
	 * @param filename - ACD file.
	 */
	public ParsedVVSADL analyseAcd() {

		// Internal application representation
		ParsedVVSADL internalRepr = null;
		
		// Define variable map
		LinkedHashMap<String, String> variableMap = new LinkedHashMap<String, String>();

		try {
			// Read the file
			File inputFile = new File(filename);
			BufferedInputStream inputStream = new BufferedInputStream(new FileInputStream(inputFile));
			final byte [] bytes = new byte[(int) inputFile.length()];
			inputStream.read(bytes);

			// ACD file analysis
			ParseAcd parser = new ParseAcd(new String(bytes), false);
			Integer numFields = parser.getNumofFields();

			// Find dependent fields (Jemboss API documentation recommends doing it)
			parser.isDependents(null, 0, numFields);

			// Read the application header
            Integer numArgs = parser.getNumofParams(0);
            for (int i = 0; i < numArgs; i++) {
                String key = parser.getParameterAttribute(0, i);
                String val = parser.getParamValueStr(0, i);
                
                // Get the key that starts with this prefix
                key = (String) getKeyByPrefix(appAttrs, key);
                if (key != null) {
                    // Store the value
                    appAttrs.put(key, val);
                }
            }
            
            // Initialize our internal representation object with header data
            internalRepr = new ParsedVVSADL(appAttrs.get("application"),
            								appAttrs.get("groups"),
            								appAttrs.get("documentation"));
            
            // Read the sections along with parameters        
            for (int j = 0; j < numFields; j++) {
                // Determine what field it is: section, endsection, parameter etc.
                String fieldType = parser.getParameterAttribute(j, 0);
                String fieldName = parser.getParamValueStr(j, 0);
                
                // Whether we should show in GUI
                Boolean showInGUI = false;

                if (fieldType.equals("section")) {
                    // A new section starts
                    // TODO: group the fields for later use in GUI

                } else if (!("endsection".startsWith(fieldType) || "application".startsWith(fieldType))) {
                    // A parameter description
                    Integer numAttrs = parser.getNumofParams(j);

                    // Loop through the attributes
                    // TODO move exp resolution to createParameter
                    for (int k = 1; k < numAttrs; k++) {
                        String attrName = parser.getParameterAttribute(j, k);
                        String attrValue = parser.getParamValueStr(j, k);

                        // Try to resolve dependencies where possible 
                        //     e.g. $(param.attr) and @($(varname)+1)
                        attrValue = resolveExp(attrValue, variableMap);

                        // Store each parameter in the variable map
                        variableMap.put(fieldName + "." + attrName, attrValue);
                        
                        // Decide whether we should show it in GUI
                        showInGUI = showInGUI ||
                                    (attrName.equals("parameter") && attrValue.equals("Y")) ||
                                    (attrName.equals("standard") && attrValue.equals("Y")) ||
                                    (attrName.equals("additional") && attrValue.equals("Y"));
                    }
                }
                
                // Check if we need to show it in GUI
                if (showInGUI) {
                    ParameterCreator.createAndAdd(parser, j, internalRepr);
                }
            }
            
            inputStream.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        return internalRepr;
    }

    /**
     * Resolve a given ACD expression.<br><br>
     * 
     * Example:<pre>
     * resolveExp("$(variable)")
     * resolveExp("@($(sequence.length)+1)")</pre>
     * 
     * @param exp - expression to be resolved.
     * @param map - a map of variable values.
     * @return resolved value.
     */
    public static String resolveExp(String exp, LinkedHashMap<String, String> map) {
        // Simulate the map
        // TODO: store some precalculated values in map (such as acdprotein)
        // TODO: deal with calculated values (such as sequence.length)

        // Regular expression for variable names like $(variable.name)
        RE reVar = new RE("\\$\\(([\\w.]+)\\)");

        // Regular expression for operation expressions @($(bool ? $(var) : 0))
        RE reFunc = new RE("\\@\\(([\\w.+-/:?\"]+)\\)");

        // Find a variable name match
        reVar.match(exp);
        String match = reVar.getParen(1);
        String resolvedExp = exp;
        if (match != null) {
            // Find replacement string
            String substitute = "$(" + match + ")";
            if (map.containsKey(match)) {
                substitute = map.get(match);
            }

            // Find position and change
            Integer start = reVar.getParenStart(1) - 2;
            Integer end = reVar.getParenEnd(1) + 1;
            resolvedExp = exp.substring(0, start) +
                          substitute +
                          exp.substring(end, exp.length());
        }

        // Find an operation expression match
        reFunc.match(resolvedExp);
        match = reVar.getParen(1);
        if (match != null) {
            AcdFunResolve resolver = new AcdFunResolve(resolvedExp);
            resolvedExp = resolver.getResult();
        }

        if (!(exp.equals(resolvedExp))) {
            return resolveExp(resolvedExp, map);
        } else {
            // No changes were made - stop the recursion
            return resolvedExp;
        }
    }

    /**
     * Given a HashMap find a key in that HashMap that
     * begins with some String prefix. Return the value
     * for the first key found.
     * 
     * @param map
     * @param prefix
     */
    private Object getKeyByPrefix(HashMap map, String prefix) {
        Object[] keys = map.keySet().toArray();
        for (Object key : keys) {
            if (((String) key).startsWith(prefix)) {
                return key;
            }
        }
        return null;
    }
    
    public static void main(String[] args) {
        ACDReader jt = new ACDReader("/opt/EMBOSS-6.2.0/emboss/acd/twofeat.acd");
        jt.analyseAcd();
    }
}