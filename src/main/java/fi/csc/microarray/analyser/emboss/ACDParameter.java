package fi.csc.microarray.analyser.emboss;

import java.util.Arrays;
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
    
    // ACDDescription that owns this parameter
    ACDDescription acd = null;
    
    // Parameter groups according to ACD specification
    static final Integer PARAM_GROUP_SIMPLE = 0;
    static final Integer PARAM_GROUP_INPUT = 1;
    static final Integer PARAM_GROUP_LIST = 2;
    static final Integer PARAM_GROUP_OUTPUT = 3;
    static final Integer PARAM_GROUP_GRAPHICS = 4;
    
    private String type;
    private String name;
    private String section;
    private String subsection;
    private LinkedHashMap<String, String> list;
    
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
     * This constructor binds this parameter to ACDDescription
     * to which this parameter belongs. There is some functionality
     * that expects parameter to know its parent ACDDescription.
     * 
     * @param acd
     * @param type
     * @param name
     * @param section
     * @param subsection
     */
    public ACDParameter(ACDDescription acd, String type, String name, String section,
                        String subsection) {
        this(type, name, section, subsection);
        this.acd = acd;
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
     * Reevaluate attributes for this parameter. E.g.
     * user has filled in some parameters on which some
     * attribute might depend on (e.g. $(paramname))
     * 
     * @param varMap - map with key/value pairs for ACD variables.
     */
    public void updateAttributes(LinkedHashMap<String, String> varMap) {
       for (String key : attributes.keySet()) {
           attributes.put(key, resolveExp(attributes.get(key), varMap));
       }
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
     * Check if given attribute has been set. An empty
     * attriute (e.g. default: "") is considered not set.
     * 
     * @param name
     * @param isEvaluated 9 - whether this parameter has to be
     *        without @(x) and $(x) stuff
     * @return true if this attribute is present, false otherwise.
     */
    public Boolean hasAttribute(String name, Boolean isEvaluated) {
        isEvaluated = !isEvaluated || attributeIsEvaluated(name);
        if (attributes.containsKey(name) &&
            !getAttribute(name).equals("") &&
            isEvaluated) {
            return true;
        }
        return false;
    }
    
    /**
     * Check if given attribute can be evaluated to True.
     * 
     * @param name of the attribute to check
     * @return
     */
    public Boolean attributeIsTrue(String name) {
        String attrValue = getAttribute(name);
        if (attrValue != null) {
            return attrValue.toLowerCase().equals("y") ||
                   attrValue.toLowerCase().equals("true");
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
        list = new LinkedHashMap<String, String>();

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
        list = new LinkedHashMap<String, String>();

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
     * Return a recommended filename for an output parameter. Valid
     * only if parameter is of some output or graph type.
     * 
     * @param withExtension - append file extension to filename.
     * @return the recommended name for output file.
     */
    public String getOutputFilename(Boolean withExtension) {
        
        // Add extension
        String extension = "";
        if (withExtension) {
            if (ACDParameter.detectParameterGroup(getType()) ==
                ACDParameter.PARAM_GROUP_OUTPUT) {
                extension = ".txt";
            } else if (ACDParameter.detectParameterGroup(getType()) ==
                       ACDParameter.PARAM_GROUP_GRAPHICS) {
                // Somehow applications adds this ".1" suffix for png files
                extension = ".1.png";
            }
        }
        
        // Add suffix
        String suffix = "";
        if (acd != null) {
            suffix = acd.getName();
        }
        
        return getName() + "." + suffix + extension;
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
        Boolean attrStandard = attributeIsTrue("standard");
        Boolean attrParameter = attributeIsTrue("parameter");
        if (attrStandard || attrParameter) {
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
        return attributeIsTrue("additional");
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
            
            // Hack booleans (somehow Jemboss does not understand Y/N)
            if (substitute.toLowerCase().equals("y")) {
                substitute = "true";
            } else if (substitute.toLowerCase().equals("n")) {
                substitute = "false";
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
     * Detect functional group of a parameter: simple, input,
     * selection list, output or graphics.
     * 
     * @param fieldType
     */
    public static Integer detectParameterGroup(String fieldType) {
        String typesSimple[] = {"array", "boolean", "float", "integer",
                                "range", "string", "toggle"};
        String typesInput[] = {"codon", "cpdb", "datafile", "directoty", "dirlist",
                               "discretestates", "distances", "features", "filelist",
                               "frequencies", "infile", "matrix", "matrixf", "pattern",
                               "properties", "regexp", "scop", "sequence", "seqall", "seqset",
                               "seqsetall", "seqsetall"};
        String typesList[] = {"list", "selection"};
        String typesOutput[] = {"align", "featout", "outcodon", "outcpdb", "outdata",
                                "outdir", "outdiscrete", "outdistance", "outfile", "outfileall",
                                "outfreq", "outmatrix", "outmatrixf", "outproperties", "outscop",
                                "outtree", "report", "seqout", "seqoutall", "seqoutset"};
        String typesGraphics[] = {"graph", "xygraph"};
        
        if (Arrays.asList(typesSimple).contains(fieldType)) {
            return PARAM_GROUP_SIMPLE;
        } else if (Arrays.asList(typesInput).contains(fieldType)) {
            return PARAM_GROUP_INPUT;
        } else if (Arrays.asList(typesList).contains(fieldType)) {
            return PARAM_GROUP_LIST;
        } else if (Arrays.asList(typesOutput).contains(fieldType)) {
            return PARAM_GROUP_OUTPUT;
        } else if (Arrays.asList(typesGraphics).contains(fieldType)) {
            return PARAM_GROUP_GRAPHICS;
        } else {
            return -1;
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
    public String normalize(String value) {
        String type = getType();
        if (type.equals("boolean")) {
            value = value.toLowerCase();
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
                
                if (hasAttribute("minimum", true)) {
                    Integer minVal = Integer.parseInt(getAttribute("minimum"));
                    accepts = accepts && (intVal >= minVal);
                }
                
                if (hasAttribute("maximum", true)) {
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
                
                if (hasAttribute("minimum", true)) {
                    Float minVal = Float.parseFloat(getAttribute("minimum"));
                    accepts = accepts && (floatVal >= minVal);
                }
                
                if (hasAttribute("maximum", true)) {
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
