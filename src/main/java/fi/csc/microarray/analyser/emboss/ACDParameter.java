package fi.csc.microarray.analyser.emboss;

import java.util.HashMap;
import java.util.LinkedHashMap;

import org.apache.regexp.RE;
import org.emboss.jemboss.parser.AcdFunResolve;

/**
 * Represents a single ACD parameter (simple, input, output etc.)
 * 
 * @author naktinis
 *
 */

public class ACDParameter {
    
    private String type;
    private String name;
    private String section;
    private String subsection;
    
    private HashMap<String, String> attributes = new HashMap<String, String>();
    
    public ACDParameter(String type, String name, String section,
                        String subsection) {
        this.setName(name);
        this.setType(type);
        this.setSection(section);
        this.setSubsection(subsection);
    }
    
    /**
     * Set parameter name.
     * 
     * @param name
     */
    public void setName(String name) {
        this.name = name;
    }
    
    /**
     * Set parameter type.
     * 
     * @param type
     */
    public void setType(String type) {
        this.type = type;
    }
    
    /**
     * Get parameter name that should be used in the
     * command line.
     * 
     * @return
     */
    public String getName() {
        return name;
    }
    
    /**
     * Get parameter type (string, integer, sequence etc.)
     * 
     * @return
     */
    public String getType() {
        return type;
    }
       
    /**
     * Define an attribute for this parameter.
     * 
     * @param name
     * @param value
     */
    public void addAttribute(String name, String value) {
        attributes.put(name, value);
    }
    
    /**
     * Find an attribute of this parameter.
     * 
     * @param name
     * @return attribute value or null if not found.
     */
    public String getAttribute(String name) {
        if (attributes.containsKey(name)) {
            return attributes.get(name);
        }
        return null;
    }
    
    /**
     * Set a section for parameter.
     * 
     * @param name
     */
    public void setSection(String name) {
        this.section = name;
    }
    
    /**
     * Get a section for parameter.
     * 
     * @return
     */
    public String getSection() {
        return section;
    }
    
    /**
     * Set a subsection for parameter.
     * 
     * @param name
     */
    public void setSubsection(String name) {
        this.subsection = name;
    }
    
    /**
     * Get a subsection for parameter.
     * 
     * @return
     */
    public String getSubsection() {
        return subsection;
    }
    
    // TODO
    
    /**
     * Determine if this parameter is required.
     * 
     * @return
     */
    public Boolean isRequired() {
        return false;
    }
    
    /**
     * Determine if this parameter is optional, but should
     * be provided for a user.
     * 
     * @return
     */
    public Boolean isAdditional() {
        return false;
    }
    
    /**
     * Determine if this parameter is advanced and should be normally
     * not shown to a user.
     * 
     * @return
     */
    public Boolean isAdvanced() {
        return false;
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
        // TODO: store some precalculated values in map (such as acdprotein) -> server side
        // TODO: deal with calculated values (such as sequence.length) -> server side

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
}
