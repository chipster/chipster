package fi.csc.microarray.analyser.emboss;

import java.util.HashMap;
import java.util.LinkedHashMap;

import org.apache.regexp.RE;
import org.emboss.jemboss.parser.AcdFunResolve;
import org.emboss.jemboss.parser.ParseAcd;

/**
 * Represents a single ACD parameter (simple, input, output, list etc.)
 * 
 * @author naktinis
 *
 */

public class ACDParameter {
    
    private String type;
    private String name;
    private String section;
    private String subsection;
    private HashMap<String, String> list;
    
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
     * Define or update an attribute for this parameter.
     * 
     * @param name
     * @param value
     */
    public void setAttribute(String name, String value) {
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
     * Define keys and values for a list parameter. Valid only
     * if paramter is of type "list" or "selection".
     * 
     * @param parser - Jemboss ParseAcd object.
     * @param index - parameter index in parser object.
     */
    public void setList(ParseAcd parser, Integer index) {
        list = new HashMap<String, String>();
        String[] values = parser.getList(index);
        
        for (int i = 0; i < values.length; i++) {
            list.put(parser.getListLabel(index, i), values[i]);
        }
    }
    
    /**
     * Get the list HashMap for this Paramter. Valid only if
     * parameter is of type "list" of "selection".
     * 
     * @return list object.
     */
    public HashMap<String, String> getList() {
        return list;
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
    
    /**
     * Determine if this parameter is required.
     * 
     * @return
     */
    public Boolean isRequired() {
        String attrStandard = getAttribute("standard");
        String attrParameter = getAttribute("parameter");
        if ((attrStandard != null && attrStandard.equals("Y")) ||
            (attrParameter != null && attrParameter.equals("Y"))) {
            return true;
        }
        return false;
    }
    
    /**
     * Determine if this parameter is optional, but should
     * be provided for a user.
     * 
     * @return
     */
    public Boolean isAdditional() {
        String attrAdditional = getAttribute("additional");
        return attrAdditional != null && attrAdditional.equals("Y");
    }
    
    /**
     * Determine if this parameter is advanced and should be normally
     * not shown to a user.
     * 
     * @return
     */
    public Boolean isAdvanced() {
        return getAttribute("standard") == null &&
               getAttribute("parameter") == null &&
               getAttribute("additional") == null;
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
