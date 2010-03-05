package fi.csc.microarray.analyser.emboss;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import org.emboss.jemboss.parser.ParseAcd;

/**
 * Represents a single ACD file.
 * 
 * @author naktinis
 * 
 */

public class ACDDescription {
    
    private String name;
    private String description;
    private LinkedList<String> groups = new LinkedList<String>(); 
    private LinkedList<ACDParameter> parameters = new LinkedList<ACDParameter>();
    private LinkedHashMap<String, String> variableMap = new LinkedHashMap<String, String>();
    
    /**
     * Get the name of this ACD.
     */
    public String getName() {
        return name;
    }
    
    /**
     * Get short informational message of this ACD.
     */
    public String getDescription() {
        return description;
    }
    
    /**
     * Get a linked list of groups this application belongs to.
     */
    public LinkedList<String> getGroups() {
        return groups;
    }
    
    /**
     * Add a parameter.
     * 
     * @param param - parameter to be added.
     */
    public void addParameter(ACDParameter param) {
        parameters.add(param);
    }
    
    /**
     * Find a parameter with given name.
     * 
     * @param name - name of the parameter to be searched for.
     * @return the parameter or null if not found.
     */
    public ACDParameter getParameter(String name) {
        for (ACDParameter parameter : parameters) {
            if (parameter.getName().equals(name)) {
                return parameter;
            }
        }
        return null;
    }
    
    /**
     * Find output parameters.
     * 
     * @return a list containing all output parameters.
     */
    public List<ACDParameter> getOutputParameters() {
        List<ACDParameter> params = new LinkedList<ACDParameter>();
        for (ACDParameter parameter : parameters) {
            if ((ACDParameter.detectParameterGroup(parameter.getType()) ==
                 ACDParameter.PARAM_GROUP_OUTPUT) &&
                parameter.isRequired()) {
                params.add(parameter);
            }
        }
        return params;
    }
    
    /**
     * Find graphics parameters.
     * 
     * @return a list containing all graphics parameters.
     */
    public List<ACDParameter> getGraphicsParameters() {
        List<ACDParameter> params = new LinkedList<ACDParameter>();
        for (ACDParameter parameter : parameters) {
            if ((ACDParameter.detectParameterGroup(parameter.getType()) ==
                 ACDParameter.PARAM_GROUP_GRAPHICS) &&
                parameter.isRequired()) {
                params.add(parameter);
            }
        }
        return params;
    }
    
    /**
     * Find all parameters in a given section or in a
     * section's subsection. You can omit the subsection.
     * Additionally, you can set <em>recursive</em> to true if you
     * want to find parmeters that are in some section as well
     * as in all its subsections.
     * 
     * @param section - name of the section.
     * @param subsection - name of the subsection or null if
     *                     you don't care about the subsection.
     * @param recursive - should we look into subsections all (valid only
     *                    if subsection is not provided).
     * @return a list of parameters that fulfil given query.
     */
    public LinkedList<ACDParameter> getParameters(String section,
                                                  String subsection,
                                                  Boolean recursive) {
        LinkedList<ACDParameter> list = new LinkedList<ACDParameter>();
        for (ACDParameter param : parameters) {
            if (param.getSection().equals(section) &&
                (subsection == null || subsection.equals(param.getSubsection())) &&
                (subsection != null || recursive || param.getSubsection() == null)) {
                list.add(param);
            }
        }
        return list;
    }
    
    /**
     * Return all parameters in this ACD.
     * 
     * @return a linked list containing all parameters.
     */
    public LinkedList<ACDParameter> getParameters() {
        return parameters;
    }
    
    /**
     * Empty constructor.
     */
    public ACDDescription() {
        
    }
    
    /**
     * Read an ACD file and store it in this object.
     * 
     * @param file - ACD file to be read.
     */
    public ACDDescription(File file) {
        
        // Application-level attributes that we want to extract
        HashMap<String, String> appAttrs = new HashMap<String, String>();
        appAttrs.put("application", "");
        appAttrs.put("documentation", "");
        appAttrs.put("groups", "");        

        try {
            // Read the file
            BufferedInputStream inputStream = new BufferedInputStream(new FileInputStream(file));
            final byte [] bytes = new byte[(int) file.length()];
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
            
            // Write some general data to the ACD object
            this.name = appAttrs.get("application");
            this.description = appAttrs.get("documentation");
            
            // Store groups
            String[] groups = appAttrs.get("groups").split(",");
            for (String group : groups) {
                this.groups.add(group.trim());
            }
            
            // Read the sections along with parameters
            String currentSection = null;
            String currentSubsection = null;
            for (int j = 1; j < numFields; j++) {
                // Determine what field it is: section, endsection, parameter etc.
                String fieldType = parser.getParameterAttribute(j, 0);
                String fieldName = parser.getParamValueStr(j, 0);

                if ("section".startsWith(fieldType)) {
                    // A new section starts
                    if (currentSection != null) {
                        currentSubsection = fieldName;
                    } else {
                        currentSection = fieldName;
                    }
                    
                } else if ("endsection".startsWith(fieldType)) {
                    // A section ends
                    if (currentSubsection == null) {
                        currentSection = null;
                    }
                    currentSubsection = null;
                    
                } else if (!("application".startsWith(fieldType))) {
                    // Initialize the parameter
                    ACDParameter param = new ACDParameter(this, fieldType, fieldName,
                            currentSection, currentSubsection);
                    
                    // A parameter description
                    Integer numAttrs = parser.getNumofParams(j);
                    
                    // Handle the lists exceptionally
                    if (fieldType.equals("list") || fieldType.equals("selection")) {
                        param.setList(parser, j);
                    }

                    // Loop through the attributes
                    for (int k = 1; k < numAttrs; k++) {
                        String attrName = parser.getParameterAttribute(j, k);
                        String attrValue = parser.getParamValueStr(j, k);
                        
                        // Try to resolve dependencies where possible 
                        //     e.g. $(param.attr) and @($(varname)+1)
                        attrValue = ACDParameter.resolveExp(attrValue, variableMap);
                        
                        // Store the attribute in the parameter
                        param.setAttribute(attrName, attrValue);

                        // Store the attribute in the variable map
                        variableMap.put(fieldName + "." + attrName, attrValue);
                    }
                    
                    // Add the parameter
                    this.addParameter(param);
                }
            }
            
            inputStream.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Update current known ACD variable values. E.g. user
     * has filled the parameters and sent them back to the
     * server.
     */
    public void updateVariables(HashMap<String, String> varMap) {
        // Update current variable map
        variableMap.putAll(varMap);
        
        // Recalculate parameter attributes
        for (ACDParameter param : getParameters()) {
            param.updateAttributes(variableMap);
        }
    }

    /**
     * Utility method. Given a HashMap find a key in that
     * HashMap that begins with some String prefix. Return
     * the value for the first key found.
     * 
     * @param map
     * @param prefix
     */
    public static Object getKeyByPrefix(HashMap<String, String> map, String prefix) {
        Object[] keys = map.keySet().toArray();
        for (Object key : keys) {
            if (((String) key).startsWith(prefix)) {
                return key;
            }
        }
        return null;
    }
}
