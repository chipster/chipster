package fi.csc.microarray.analyser.emboss;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
    
    private final HashMap<String, ACDValidator> validators =
        new HashMap<String, ACDValidator>();
    
    public ACDParameter(String type, String name, String section,
                        String subsection) {
        this.setName(name);
        this.setType(type);
        this.setSection(section);
        this.setSubsection(subsection);
        
        validators.put("float", new FloatValidator());
        validators.put("integer", new IntegerValidator());
        validators.put("array", new ArrayValidator());
        validators.put("boolean", new BooleanValidator());
        validators.put("list", new ListValidator());
        validators.put("selection", new ListValidator());
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
        this.type = type.toLowerCase();
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
     * Check if given attribute has been set. An empty
     * attriute (e.g. default: "") is considered not set.
     * 
     * @param name
     * @return true if this attribute is present, false otherwise.
     */
    public Boolean hasAttribute(String name) {
        if (attributes.containsKey(name) && !getAttribute(name).equals("")) {
            return true;
        }
        return false;
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

        // Check if it is "list" or "selection"
        if (getType().equals("list")) {
            String[] titles = parser.getList(index);
            
            for (int i = 0; i < titles.length; i++) {
                list.put(parser.getListLabel(index, i), titles[i]);
            }
        } else {
            String[] titles = parser.getSelect(index);

            for (Integer i = 0; i < parser.getSelect(index).length; i++) {
                list.put(i.toString(), titles[i]);
            }
        }
    }
    
    /**
     * Define keys and values for a list parameter. Valid only
     * if paramter is of type "list" or "selection".
     * 
     * @param titles - array containing titles.
     * @param values - array containing values for corresponding titles;
     * ignored for "selection" type.
     */
    public void setList(String[] titles, String[] values) {
        list = new HashMap<String, String>();

        // Check if it is "list" or "selection"
        if (getType().equals("list")) {           
            for (int i = 0; i < titles.length; i++) {
                list.put(values[i], titles[i]);
            }
        } else {
            for (Integer i = 0; i < titles.length; i++) {
                list.put(i.toString(), titles[i]);
            }
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
     * Determine if some attribute does not require additional
     * evaluation. If parameter has a value such as "$(varname)"
     * or "@($(varname)+1)" it is considered unevaluated.
     * 
     * @param attrName - attribute to be checked.
     * @return true if attribute value does not require further
     * evaluation.
     */
    public Boolean attributeIsEvaluated(String attrName) {
        String attrValue = getAttribute(attrName);
        return !(attrValue == null || attrValue.contains("$") || attrValue.contains("@"));
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
    
    /**
     * Normalize given value according to this parameter.
     * For example, if the type is "boolean", we should
     * convert "true" to "Y" (i.e. we convert the value
     * according to ACD spec.)
     * 
     * This method should be used to convert values,
     * entered in GUI to values that can be passed to
     * actual EMBOSS application in command line.
     * 
     * @param value - ACD-unaware value.
     */
    // TODO check out what we get for ENUM, arrays etc.
    public String normalize(String value) {
        String type = getType();
        value = value.toLowerCase();
        if (type.equals("boolean")) {
            if (value.equals("yes") || value.equals("true") || value.equals("1") ||
                value.equals("y")) {
                return "Y";
            } else if (value.equals("no") || value.equals("true") || value.equals("0") ||
                       value.equals("n")) {
                return "N";
            }
        } else if (type.equals("selection") || type.equals("list")) {
            String defaultDelim = ",";
            String[] choices = value.split(defaultDelim);
            String attrDelim = hasAttribute("delimiter") ? getAttribute("delimiter") : ";";
            
            // Reconnect values using different delimiter
            if (choices.length > 1 && attrDelim != defaultDelim) {
                String normalValue = "";
                for (String choice : choices) {
                    normalValue = normalValue.concat(attrDelim).concat(choice);
                }
                return normalValue.substring(1);
            }
        }
        return value;
    }
    
    /**
     * Check if a given value can be passed as this
     * parameter.
     * 
     * @param value
     */
    public boolean validate(String value) {
        ACDValidator validator;
        if (validators.containsKey(getType())) {
            validator = validators.get(getType());
        } else {
            validator = new VoidValidator();
        }

        return validator.accepts(normalize(value));
    }
    
    // TODO: validate ranges, input, output files
    
    /**
     * Used for data types which are not (not yet) validated.
     */
    class VoidValidator extends ACDValidator {       
        public boolean accepts(String value) {
            return true;
        }
    }
    
    class BooleanValidator extends ACDValidator {       
        public boolean accepts(String value) {
            if (value == "Y" || value == "N") {
                return true;
            }
            return false;
        }
    }
    
    class ArrayValidator extends ACDValidator {       
        public boolean accepts(String value) {
            Pattern re = Pattern.compile("^(\\d+(\\.\\d+)?[ ,])*\\d+(\\.\\d+)?$");
            Matcher m = re.matcher(value);
            if (!m.find()) {
                return false;
            }
            return true;
        }
    }
    
    class IntegerValidator extends ACDValidator {       
        public boolean accepts(String value) {
            Boolean accepts = true;
            try {
                Integer intVal = Integer.parseInt(value);
                
                if (hasAttribute("minimum")) {
                    Integer minVal = Integer.parseInt(getAttribute("minimum"));
                    accepts = accepts && (intVal >= minVal);
                }
                
                if (hasAttribute("maximum")) {
                    Integer maxVal = Integer.parseInt(getAttribute("maximum"));
                    accepts = accepts && (intVal <= maxVal);
                }
                return accepts;
            }
            catch(NumberFormatException nfe) {
                return false;
            }
        }
    }
    
    class FloatValidator extends ACDValidator {       
        public boolean accepts(String value) {
            Boolean accepts = true;
            try {
                Float floatVal = Float.parseFloat(value);
                
                if (hasAttribute("minimum")) {
                    Float minVal = Float.parseFloat(getAttribute("minimum"));
                    accepts = accepts && (floatVal >= minVal);
                }
                
                if (hasAttribute("maximum")) {
                    Float maxVal = Float.parseFloat(getAttribute("maximum"));
                    accepts = accepts && (floatVal <= maxVal);
                }                
                return accepts;
            }
            catch(NumberFormatException nfe) {
                return false;
            }
        }
    }
    
    class ListValidator extends ACDValidator {       
        public boolean accepts(String value) {
            Boolean accepts = true;
            
            HashMap<String, String> acceptedList = getList();
            String attrDelim = hasAttribute("delimiter") ? getAttribute("delimiter") : ";";
            String[] choices = value.split(attrDelim);
            
            if (hasAttribute("minimum")) {
                Integer minVal = Integer.parseInt(getAttribute("minimum"));
                accepts = accepts && (choices.length >= minVal);
            }
            
            if (hasAttribute("maximum")) {
                Integer maxVal = Integer.parseInt(getAttribute("maximum"));
                accepts = accepts && (choices.length <= maxVal);
            }
            
            for (String choice : choices) {
                accepts = accepts && acceptedList.containsKey(choice);
            }

            return accepts;
        }
    }
    
    private abstract class ACDValidator {       
        abstract public boolean accepts(String value);
    }
}
