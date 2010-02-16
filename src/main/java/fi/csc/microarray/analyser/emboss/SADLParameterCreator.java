package fi.csc.microarray.analyser.emboss;

import java.util.Arrays;
import java.util.HashMap;

import org.emboss.jemboss.parser.ParseAcd;

import fi.csc.microarray.description.GenericInputTypes;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.description.SADLSyntax.ParameterType;

/**
 * Create SADL parameters using information from ACD parser.
 * 
 * @author naktinis
 * 
 */
public class SADLParameterCreator {
    
    private static final Integer PARAM_GROUP_SIMPLE = 0;
    private static final Integer PARAM_GROUP_INPUT = 1;
    private static final Integer PARAM_GROUP_LIST = 2;
    private static final Integer PARAM_GROUP_OUTPUT = 3;
    private static final Integer PARAM_GROUP_GRAPHICS = 4;
    
    /**
     * Create a SADL parameter (Parameter, Input or Output) and
     * add it to a given object.
     * 
     * @param parser - ACD Parser object (to read the parameter).
     * @param index - parameter index in the ACD Parser (to read the parameter).
     * @param internalRepr - the object to add this parameter to.
     */
    public static void createAndAdd(ParseAcd parser, Integer index, SADLDescription internalRepr) {
        // Try to create a parameter (simple types, list types)
        Parameter param = createParameter(parser, index);
        if (param != null) {
            internalRepr.addParameter(param);
            return;
        }
        
        // Try to create an input
        Input input = createInput(parser, index);
        if (input != null) {
            internalRepr.addInput(input);
            return;
        }
        
        // Try to create an output
        String output = createOutput(parser, index);
        if (output != null) {
            internalRepr.addOutput(output);
        }
    }
    
    /**
     * Create a parameter and add it to a given SADL object.
     * Creates simple parameters, such as integers, strings etc.
     * 
     * @param parser - parser object.
     * @param index - index of a field to be parsed.
     * @return vvsadl parameter object or null.
     */
    public static Parameter createParameter(ParseAcd parser, Integer index) {
        String fieldType = parser.getParameterAttribute(index, 0);
        String fieldName = parser.getParamValueStr(index, 0);
        
        // Detect the parameter functional group
        Integer type = detectParameterGroup(fieldType.toLowerCase());
        
        // Map simple ACD parameters to VVSADL parameters
        HashMap<String, ParameterType> typeMap = new HashMap<String, ParameterType>();
        typeMap.put("array", ParameterType.ENUM);
        typeMap.put("float", ParameterType.DECIMAL);
        typeMap.put("integer", ParameterType.INTEGER);
        typeMap.put("string", ParameterType.STRING);
        typeMap.put("range", ParameterType.STRING);
        
        // Read common attributes
        String fieldDefault = parser.getDefaultParamValueStr(index);
        // TODO: help attribute; comment attribute
        String fieldInfo = parser.getInfoParamValue(index);
        
        if (fieldType == "boolean" || fieldType == "toggle") {
            // Boolean types need some special handling
            String[] fieldOptions = {"Yes", "No"};
            return new Parameter(fieldName, typeMap.get(fieldType), fieldOptions,
                    null, null, fieldDefault, fieldInfo);
        } else if (type == PARAM_GROUP_SIMPLE) {
            String fieldMin = parser.getMinParam(index);
            String fieldMax = parser.getMaxParam(index);
            return new Parameter(fieldName, typeMap.get(fieldType), null,
                                 fieldMin, fieldMax, fieldDefault, fieldInfo);
        } else if (type == PARAM_GROUP_LIST) {
            String[] fieldOptions = parser.getList(index);
            // TODO: lists with labels
            // System.out.println(parser.getListLabel(index, 1));
            return new Parameter(fieldName, typeMap.get(fieldType), fieldOptions,
                                 null, null, fieldDefault, fieldInfo);
        } else {
            return null;
        }
    }
    
    /**
     * Create a SADL input and add it to a given object.
     * Creates objects representing input files.
     * 
     * @param parser - parser object.
     * @param index - index of a field to be parsed.
     * @return vvsadl parameter object or null.
     */
    public static Input createInput(ParseAcd parser, Integer index) { 
        String fieldType = parser.getParameterAttribute(index, 0);
        String fieldName = parser.getParamValueStr(index, 0);
        
        // Detect the parameter functional group
        Integer type = detectParameterGroup(fieldType.toLowerCase());
        
        // TODO: help attribute; comment attribute
        // String fieldInfo = parser.getInfoParamValue(index);
        
        if (type == PARAM_GROUP_INPUT) {
            return Input.createInput(GenericInputTypes.GENERIC, fieldName);
        } else {
            return null;
        }
    }
    
    /**
     * Create objects representing output files (currently
     * returns a string).
     * 
     * @param parser - parser object.
     * @param index - index of a field to be parsed.
     * @return vvsadl parameter object or null.
     */
    public static String createOutput(ParseAcd parser, Integer index) { 
        String fieldType = parser.getParameterAttribute(index, 0);
        String fieldName = parser.getParamValueStr(index, 0);
        
        // Detect the parameter functional group
        Integer type = detectParameterGroup(fieldType.toLowerCase());
        
        // TODO: help attribute; comment attribute
        // String fieldInfo = parser.getInfoParamValue(index);
        
        if (type == PARAM_GROUP_OUTPUT) {
            return fieldName;
        } else {
            return null;
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
}
