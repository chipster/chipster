package fi.csc.microarray.description;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.description.SADLDescription.Name;
import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLDescription.Output;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.description.SADLSyntax.InputType;
import fi.csc.microarray.description.SADLSyntax.ParameterType;
import fi.csc.microarray.description.vvsadl.CompatibilityVVSADLParser;
import fi.csc.microarray.exception.MicroarrayException;



/**
 * <p>Parses SADL analysis descriptions. SADL stands for  
 * Simple Analysis Description Language. It is used to manage analysis 
 * description data inside the Chipster system. Parsing is event based (cf. SAX).</p>
 *  
 * @see fi.csc.microarray.description.SADLSyntax
 * 
 * @author Aleksi Kallio
 */
public class SADLParser {
		
	/**
	 * Logger for this class (debugging).
	 */
	private static final Logger logger = Logger.getLogger(SADLParser.class);

	
	/**
	 * Parse failure caused by illegal input data (bad SADL).
	 */
	public static class ParseException extends MicroarrayException {
		public ParseException(String msg) {
			super(msg);
		}
		
		public ParseException(String msg, String filename) {
			super(msg + (filename != null ? (" (in " + filename + ")") : ""));
		}
	}
	
	private String unitName;
	private HashMap<String, InputType> inputTypeMap = new HashMap<String, InputType>();
	
	public static String generateOperationIdentifier(String category, String name) {
		return "\"" + category + "\"/\"" + name + "\"";
	}
	
	public SADLParser() {
		this(null);
	}
	
	public SADLParser(String filename) {
		this.unitName = filename;
		addInputType(GenericInputTypes.GENERIC);
	}

	public SADLDescription parse(String sadlString) throws ParseException {
		
		// check for VVSADL compatibility mode
		if (sadlString.trim().startsWith("ANALYSIS")) {
			return new CompatibilityVVSADLParser().parse(sadlString);
		}
		
		SADLTokeniser tokens = new SADLTokeniser(sadlString, unitName);
		return parseTool(tokens);
	}
	
    public SADLDescription parse(String sadlString, String id) throws ParseException {
        
        // check for VVSADL compatibility mode
        if (sadlString.trim().startsWith("ANALYSIS")) {
            SADLDescription sadl = new CompatibilityVVSADLParser().parse(sadlString);
            sadl.setID(id);
            return sadl;
        }
        
        SADLTokeniser tokens = new SADLTokeniser(sadlString, unitName);
        return parseTool(tokens);
    }
	
	public List<SADLDescription> parseMultiple(String sadlString) throws ParseException {
		
		// check for VVSADL compatibility mode
		if (sadlString.trim().startsWith("ANALYSIS")) {
			return new CompatibilityVVSADLParser().parseMultiple(sadlString);			
		}
		
		LinkedList<SADLDescription> descriptions = new LinkedList<SADLDescription>();
		SADLTokeniser tokens = new SADLTokeniser(sadlString, unitName);
		
		// for avoiding excessive logging in case of parse error
		boolean parsingPreviousSuccessful = true;
		while (tokens.hasNext()) {
			try {
				descriptions.add(parseTool(tokens));
				parsingPreviousSuccessful = true;
			} catch (ParseException pe) {
				if (parsingPreviousSuccessful) {
					descriptions.removeLast();
					logger.error("Could not parse description. ", pe);
					
				} 
				parsingPreviousSuccessful = false;
			}
		}
	
		return descriptions;
	}


	public void addInputType(InputType type) {
		this.inputTypeMap .put(type.getName(), type);
	}
	
	/**
	 * Parsing is implemented with recursive descent algorithm (with 1 token look-a-head).
	 */
	private SADLDescription parseTool(SADLTokeniser tokens) throws ParseException {
		// read first line (analysis)
		skip(tokens, SADLSyntax.KEYWORD_TOOL);

		// read analysis stuff
		String category = tokens.next();
		skip(tokens, SADLSyntax.CATEGORY_SEPARATOR);
		Name name = parseName(tokens);		
		String comment = tokens.next();
		SADLDescription description = new SADLDescription(name, category, comment);
	
		// read possible inputs
		while (nextTokenIs(tokens, SADLSyntax.KEYWORD_INPUT)) { 
			skip(tokens, SADLSyntax.KEYWORD_INPUT);  
			Input input = parseInput(tokens, description);
			description.addInput(input);			
		}

		// read possible metainputs
		while (nextTokenIs(tokens, SADLSyntax.KEYWORD_METAINPUT)) { 
			skip(tokens, SADLSyntax.KEYWORD_METAINPUT); 
			Input input = parseInput(tokens, description);
			description.addMetaInput(input);			
		}

		// read possible outputs
		while (nextTokenIs(tokens, SADLSyntax.KEYWORD_OUTPUT)) {
			skip(tokens, SADLSyntax.KEYWORD_OUTPUT); 
			description.addOutput(parseOutput(tokens));
		}

		// read possible metaoutputs
		while (nextTokenIs(tokens, SADLSyntax.KEYWORD_METAOUTPUT)) {
			skip(tokens, SADLSyntax.KEYWORD_METAOUTPUT);
			description.addMetaOutput(parseOutput(tokens));
		}

		//	read possible parameters
		while (nextTokenIs(tokens, SADLSyntax.KEYWORD_PARAMETER)) {
			skip(tokens, SADLSyntax.KEYWORD_PARAMETER);
			Parameter parameter = parseParameter(tokens);
			description.addParameter(parameter);
		}
		
		return description;
	}

	private void skip(SADLTokeniser tokens, String token) throws ParseException {
		String next = tokens.next();
		if (!token.equals(next)) {
			throw new ParseException("expected " + token + ", not " + next, unitName);
		}
	}

	private boolean nextTokenIs(SADLTokeniser tokens, String token) {
		return tokens.hasNext() && token.equals(tokens.peek());
	}

	private Name parseName(SADLTokeniser tokens) throws ParseException {
		Name name = SADLDescription.Name.createEmptyName();
		
		String rawName = tokens.next();
		
		if (rawName.contains(SADLSyntax.NAME_SET_DESIGNATOR)) {
			name.setPrefix(rawName.substring(0, rawName.indexOf(SADLSyntax.NAME_SET_DESIGNATOR)));
			name.setPostfix(rawName.substring(rawName.indexOf(SADLSyntax.NAME_SET_DESIGNATOR) + SADLSyntax.NAME_SET_DESIGNATOR.length()));
			
		} else {
			name.setID(rawName);
		}
		
		if (tokens.hasNext() && SADLSyntax.NAME_SEPARATOR.equals(tokens.peek())) {
			skip(tokens, SADLSyntax.NAME_SEPARATOR); // read separator
			name.setDisplayName(tokens.next());
		}
		
		// check 
		if (name == null) {
			throw new ParseException("name is null");
		}
		if (name.getID() == null && (name.getPrefix() == null && name.getPostfix() == null)) {
			throw new ParseException("id, prefix and postfix are all null");
		}
		
		return name;
	}

	private Output parseOutput(SADLTokeniser tokens) throws ParseException {
		
		Output output = new Output();
		boolean isOptional = parseOptionalIfExists(tokens);
		output.setOptional(isOptional);

		output.setName(parseName(tokens));
		
		return output;
	}

	private Input parseInput(SADLTokeniser tokens, SADLDescription description) throws ParseException {
		
		Input input = new Input();
		boolean isOptional = parseOptionalIfExists(tokens);
		input.setOptional(isOptional);
		
		input.setName(parseName(tokens));
		skip(tokens, SADLSyntax.KEYWORD_TYPE);  
		String inputTypeName = tokens.next();
		InputType inputType = inputTypeMap.get(inputTypeName);
		if (inputType == null) {
			throw new ParseException("Invalid input type: " + inputTypeName, description.getName().getID());
		}
		input.setType(inputType);
				
		return input;
	}
	
	private boolean parseOptionalIfExists(SADLTokeniser tokens) throws ParseException {
		
		if (nextTokenIs(tokens, SADLSyntax.KEYWORD_OPTIONAL)) {
			skip(tokens, SADLSyntax.KEYWORD_OPTIONAL);
			return true;
		
		} else {
			return false;
		}
	}

	private Parameter parseParameter(SADLTokeniser tokens) throws ParseException {
		
		boolean isOptional = parseOptionalIfExists(tokens);
		
		Name name = parseName(tokens);
		
		skip(tokens, SADLSyntax.KEYWORD_TYPE);
		ParameterType type = null;
		Name[] options = null;
		
		if (nextTokenIs(tokens, SADLSyntax.ENUM_OPEN)) {
			options = parseEnumType(tokens);
			type = ParameterType.ENUM;
			
		} else {
			type = ParameterType.valueOf(tokens.next());
		}
		
		String from = null;
		String to = null;
		String[] defaultValues = new String[0]; 

		if (nextTokenIs(tokens, SADLSyntax.KEYWORD_FROM)) {
			skip(tokens, SADLSyntax.KEYWORD_FROM); 
			from = tokens.next(); 
		}

		if (nextTokenIs(tokens, SADLSyntax.KEYWORD_TO)) {
			skip(tokens, SADLSyntax.KEYWORD_TO); 
			to = tokens.next();
		}

		if (nextTokenIs(tokens, SADLSyntax.KEYWORD_DEFAULT)) {
			skip(tokens, SADLSyntax.KEYWORD_DEFAULT);
			defaultValues = parseDefaultValues(tokens);
		}

		String comment = tokens.next();
		
		Parameter parameter = new Parameter(name, type, options, from, to, defaultValues, comment);
		parameter.setOptional(isOptional);
		
		return parameter;
	}

	private String[] parseDefaultValues(SADLTokeniser tokens) throws ParseException {
		
		LinkedList<String> list = new LinkedList<String>();		
		while (true) {
			list.add(tokens.next());
			if (nextTokenIs(tokens, SADLSyntax.LIST_SEPARATOR)) {
				skip(tokens, SADLSyntax.LIST_SEPARATOR);
			} else {
				break;	
			}			
		}
		
		return list.toArray(new String[0]);
	}


	private Name[] parseEnumType(SADLTokeniser tokens) throws ParseException {
		
		skip(tokens, SADLSyntax.ENUM_OPEN);
		
		LinkedList<Name> list = new LinkedList<Name>();		
		while (true) {
			list.add(parseName(tokens));
			if (nextTokenIs(tokens, SADLSyntax.ENUM_CLOSE)) {
				break;	
			} else {
				skip(tokens, SADLSyntax.LIST_SEPARATOR);				
			}			
		}
		
		skip(tokens, SADLSyntax.ENUM_CLOSE);

		return list.toArray(new Name[0]);
	}
	
}