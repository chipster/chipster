package fi.csc.microarray.description;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLDescription.AnnotatedName;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.description.SADLSyntax.InputType;
import fi.csc.microarray.description.SADLSyntax.ParameterType;
import fi.csc.microarray.exception.MicroarrayException;



/**
 * <p>Parses VVSADL analysis descriptions. VVSADL stands for Very Very 
 * Simple Analysis Description Language. It is used to manage analysis 
 * description data inside the Chipster system. Parsing is event based (cf. SAX).</p>
 *  
 * @see fi.csc.microarray.description.SADLSyntax
 * 
 * @author Aleksi Kallio
 */
public class SADLParser {

	
	private static final String INPUT_SET_DESIGNATOR = "{...}";
	private static final String NAME_SEPARATOR = ":";
	/**
	 * Logger for this class (debugging).
	 */
	private static final Logger logger = Logger.getLogger(SADLParser.class);

	
	/**
	 * Parse failure caused by illegal input data (bad VVSADL).
	 */
	public static class ParseException extends MicroarrayException {
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

	public SADLDescription parse(String vvsadlString) throws ParseException {
		SADLTokeniser tokens = new SADLTokeniser(vvsadlString, unitName);
		return parseAnalysis(tokens);
	}

	public List<SADLDescription> parseMultiple(String vvsadlString) throws ParseException {
		LinkedList<SADLDescription> descriptions = new LinkedList<SADLDescription>();
		SADLTokeniser tokens = new SADLTokeniser(vvsadlString, unitName);
		
		// for avoiding excessive logging in case of parse error
		boolean parsingPreviousSuccessful = true;
		while (tokens.hasNext()) {
			try {
				descriptions.add(parseAnalysis(tokens));
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
	private SADLDescription parseAnalysis(SADLTokeniser tokens) throws ParseException {
		// read first line (analysis)
		skip(tokens, "TOOL");

		// read analysis stuff
		AnnotatedName name = parseName(tokens);
		String comment = tokens.next();
		SADLDescription description = new SADLDescription(name, comment);
	
		// read possible inputs
		while (nextTokenIs(tokens, "INPUT")) { 
			skip(tokens, "INPUT");  
			Input input = parseInput(tokens, description);
			description.addInput(input);			
		}

		// read possible metainputs
		while (nextTokenIs(tokens, "METAINPUT")) { 
			skip(tokens, "METAINPUT"); 
			Input input = parseInput(tokens, description);
			description.addMetaInput(input);			
		}

		// read possible outputs
		while (nextTokenIs(tokens, "OUTPUT")) {
			skip(tokens, "OUTPUT"); 
			String output = parseOutput(tokens);
			description.addOutput(output);
		}

		// read possible metaoutputs
		while (nextTokenIs(tokens, "METAOUTPUT")) {
			skip(tokens, "METAOUTPUT");
			String output = parseOutput(tokens);
			description.addMetaOutput(output);
		}

		//	read possible parameters
		while (nextTokenIs(tokens, "PARAMETER")) {
			skip(tokens, "PARAMETER");
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

	private AnnotatedName parseName(SADLTokeniser tokens) throws ParseException {
		AnnotatedName name = SADLDescription.AnnotatedName.createEmptyName();
		
		String rawName = tokens.next();
		
		if (rawName.contains(INPUT_SET_DESIGNATOR)) {
			name.setPrefix(rawName.substring(0, rawName.indexOf(INPUT_SET_DESIGNATOR)));
			name.setPostfix(rawName.substring(rawName.indexOf(INPUT_SET_DESIGNATOR) + INPUT_SET_DESIGNATOR.length()));
			
		} else {
			name.setName(rawName);
		}

		
		if (NAME_SEPARATOR.equals(tokens.peek())) {
			skip(tokens, NAME_SEPARATOR); // read separator
			name.setHumanReadableName(tokens.next());
		}
		return name;
	}

	private String parseOutput(SADLTokeniser tokens) throws ParseException {
		
		if (nextTokenIs(tokens, "OPTIONAL")) {
			skip(tokens, "OPTIONAL");
		}

		AnnotatedName name = parseName(tokens);
		
		return name.getName();
	}

	private Input parseInput(SADLTokeniser tokens, SADLDescription description) throws ParseException {
		
		Input newInput = new Input();
		boolean isOptional = parseOptionalIfExists(tokens);
		newInput.setOptional(isOptional);
		
		newInput.setName(parseName(tokens));
		skip(tokens, "TYPE");  
		newInput.setType(inputTypeMap.get(tokens.next()));
				
		return newInput;
	}
	
	private boolean parseOptionalIfExists(SADLTokeniser tokens) throws ParseException {
		
		if (nextTokenIs(tokens, "OPTIONAL")) {
			skip(tokens, "OPTIONAL");
			return true;
		
		} else {
			return false;
		}
	}

	private Parameter parseParameter(SADLTokeniser tokens) throws ParseException {
		
		boolean isOptional = parseOptionalIfExists(tokens);
		
		AnnotatedName name = parseName(tokens);
		
		skip(tokens, "TYPE");
		ParameterType type = null;
		String[] options = null;
		
		if (nextTokenIs(tokens, "[")) {
			options = parseEnumType(tokens);
			type = ParameterType.ENUM;
			
		} else {
			type = ParameterType.valueOf(tokens.next());
		}
		
		String from = null;
		String to = null;
		String defaultValue = null; 

		if (nextTokenIs(tokens, "FROM")) {
			skip(tokens, "FROM"); 
			from = tokens.next(); 
		}

		if (nextTokenIs(tokens, "TO")) {
			skip(tokens, "TO"); 
			to = tokens.next();
		}

		if (nextTokenIs(tokens, "DEFAULT")) {
			skip(tokens, "DEFAULT"); 
			defaultValue = tokens.next();
		}

		String comment = tokens.next();
		
		Parameter parameter = new Parameter(name, type, options, from, to, defaultValue, comment);
		
		return parameter;
	}


	private String[] parseEnumType(SADLTokeniser tokens) throws ParseException {
		
		skip(tokens, "[");
		
		LinkedList<String> list = new LinkedList<String>();		
		while (true) {
			list.add(parseName(tokens).getName());
			if (nextTokenIs(tokens, "]")) {
				break;	
			} else {
				skip(tokens, ",");				
			}			
		}
		
		skip(tokens, "]");

		return list.toArray(new String[0]);
	}
	
}