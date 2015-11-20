package fi.csc.microarray.description;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLDescription.Name;
import fi.csc.microarray.description.SADLDescription.Output;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.description.SADLSyntax.InputType;
import fi.csc.microarray.description.SADLSyntax.ParameterType;
import fi.csc.microarray.description.SADLTokeniser.TokenType;
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
	@SuppressWarnings("serial")
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

		SADLTokeniser tokens = new SADLTokeniser(sadlString, unitName);
		return parseTool(tokens);
	}
	
	public List<SADLDescription> parseMultiple(String sadlString) throws ParseException {
		
		LinkedList<SADLDescription> descriptions = new LinkedList<SADLDescription>();
		SADLTokeniser tokens = new SADLTokeniser(sadlString, unitName);
		
		// for avoiding excessive logging in case of parse error
		boolean parsingPreviousSuccessful = true;
		while (tokens.hasNext()) {
			try {
				descriptions.add(parseTool(tokens));
				parsingPreviousSuccessful = true;
			} catch (ParseException pe) {
				pe.printStackTrace();
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
		Name name = parseName(tokens);		
		SADLDescription description = new SADLDescription(name);

		if (tokens.peekType() == TokenType.DESCRIPTION) {
			description.setDescription(tokens.next());
		}

		// read possible inputs
		while (nextTokenIs(tokens, SADLSyntax.KEYWORD_INPUT)) { 
			skip(tokens, SADLSyntax.KEYWORD_INPUT);  
			Input input = parseInput(tokens, description);
			description.addInput(input);			
		}

		// read possible outputs
		while (nextTokenIs(tokens, SADLSyntax.KEYWORD_OUTPUT)) {
			skip(tokens, SADLSyntax.KEYWORD_OUTPUT); 
			description.addOutput(parseOutput(tokens));
		}
		
		//	read possible parameters
		while (nextTokenIs(tokens, SADLSyntax.KEYWORD_PARAMETER)) {
			skip(tokens, SADLSyntax.KEYWORD_PARAMETER);
			Parameter parameter = parseParameter(tokens);
			description.addParameter(parameter);
		}

		// check that no trailing content was left behind
		if (tokens.hasNext() && !nextTokenIs(tokens, SADLSyntax.KEYWORD_TOOL)) {
			// content other then new description was left 
			throw new ParseException("unexpected content: " + tokens.next(), unitName);
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
		
		if (name.getID() == null && (name.getPrefix() == null && name.getPostfix() == null)) {
			throw new ParseException("id, prefix and postfix are all null");
		}
		
		return name;
	}

	private Output parseOutput(SADLTokeniser tokens) throws ParseException {
		
		Output output = new Output();

		boolean isMeta = parseMetaIfExists(tokens);
		output.setMeta(isMeta);

		boolean isOptional = parseOptionalIfExists(tokens);
		output.setOptional(isOptional);

		output.setName(parseName(tokens));

		if (tokens.peekType() == TokenType.DESCRIPTION) {
			output.setDescription(tokens.next());
		}

		return output;
	}

	private Input parseInput(SADLTokeniser tokens, SADLDescription description) throws ParseException {
		
		Input input = new Input();

		boolean isMeta = parseMetaIfExists(tokens);
		input.setMeta(isMeta);
		
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
		
		if (tokens.peekType() == TokenType.DESCRIPTION) {
			input.setDescription(tokens.next());
		}

		return input;
	}

	private boolean parseMetaIfExists(SADLTokeniser tokens) throws ParseException {
		
		if (nextTokenIs(tokens, SADLSyntax.KEYWORD_META)) {
			skip(tokens, SADLSyntax.KEYWORD_META);
			return true;
		
		} else {
			return false;
		}
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

		Parameter parameter = new Parameter(name, type, options, from, to, defaultValues);
		parameter.setOptional(isOptional);

		if (tokens.peekType() == TokenType.DESCRIPTION) {
			parameter.setDescription(tokens.next());
		}

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