package fi.csc.microarray.description;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.description.ParsedVVSADL.Input;
import fi.csc.microarray.description.ParsedVVSADL.Parameter;
import fi.csc.microarray.description.VVSADLSyntax.InputType;
import fi.csc.microarray.description.VVSADLSyntax.ParameterType;
import fi.csc.microarray.util.AdvancedStringTokenizer;
import fi.csc.microarray.util.Deseparator;



/**
 * <p>Parses VVSADL analysis descriptions. VVSADL stands for Very Very 
 * Simple Analysis Description Language. It is used to manage analysis 
 * description data inside the Chipster system. Parsing is event based (cf. SAX).</p>
 *  
 * @see fi.csc.microarray.description.VVSADLSyntax
 * 
 * @author Aleksi Kallio
 */
public class VVSADLParser {

	
	/**
	 * Logger for this class (debugging).
	 */
	private static final Logger logger = Logger.getLogger(VVSADLParser.class);

	
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
	
	public VVSADLParser() {
		this(null);
	}
	
	public VVSADLParser(String filename) {
		this.unitName = filename;
		addInputType(GenericInputTypes.GENERIC);
	}

	public ParsedVVSADL parse(String vvsadlString) throws ParseException {
		AdvancedStringTokenizer tokens = getTokenizer(vvsadlString);
		return parseAnalysis(tokens);
	}

	public List<ParsedVVSADL> parseMultiple(String vvsadlString) throws ParseException {
		LinkedList<ParsedVVSADL> descriptions = new LinkedList<ParsedVVSADL>();
		AdvancedStringTokenizer tokens = getTokenizer(vvsadlString);
		
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
	private ParsedVVSADL parseAnalysis(AdvancedStringTokenizer tokens) throws ParseException {
		// read first line (analysis)
		if (!"ANALYSIS".equals(tokens.next())) {
			throw new ParseException("VVSADL should start with \"ANALYSIS\", not " + tokens.current(), unitName);
		}

		// read analysis stuff
		String packageName = tokens.next();
		String name = tokens.next();
		if (this.unitName == null) {
			this.unitName = packageName + " / " + name;
		}
		String comment = readComment(tokens);
		ParsedVVSADL description = new ParsedVVSADL(name, packageName, comment);
	
		// read possible inputs
		if (tokens.hasNext() && "INPUT".equals(tokens.peek())) { 
			tokens.next(); // skip "INPUT"
			List<Input> inputs = parseInputs(tokens, description);
			description.addInputs(inputs);			
		}

		// read possible metainputs
		if (tokens.hasNext() && "METAINPUT".equals(tokens.peek())) { 
			tokens.next(); // skip "METAINPUT"
			List<Input> inputs = parseInputs(tokens, description);
			description.addMetaInputs(inputs);			
		}

		// read possible outputs
		if (tokens.hasNext() && "OUTPUT".equals(tokens.peek())) {
			tokens.next(); // skip "OUTPUT"
			List<String> outputs = parseOutputs(tokens, description);
			description.addOutputs(outputs);
		}

		// read possible metaoutputs
		if (tokens.hasNext() && "METAOUTPUT".equals(tokens.peek())) {
			tokens.next(); // skip "METAOUTPUT"
			List<String> outputs = parseOutputs(tokens, description);
			description.addMetaOutputs(outputs);
		}

		//	read possible parameters
		while (tokens.hasNext() && "PARAMETER".equals(tokens.peek())) {
			Parameter parameter = parseParameter(tokens);
			description.addParameter(parameter);
		}
		
		return description;
	}

	private List<String> parseOutputs(AdvancedStringTokenizer tokens, ParsedVVSADL description) {
		LinkedList<String> outputList = new LinkedList<String>();
		Deseparator outputs = new Deseparator(",", tokens, 1);
		for (String[] output : outputs) {
			outputList.add(output[0]);
		}
		return outputList;
	}

	private List<Input> parseInputs(AdvancedStringTokenizer tokens, ParsedVVSADL description) {
		LinkedList<Input> inputList = new LinkedList<Input>();
		Deseparator inputs = new Deseparator(",", tokens, 2);
		for (String[] input : inputs) {
			Input newInput;
			if (input[1].contains("[...]")) {
				String filePattern = input[1];
				String prefix = filePattern.substring(0, filePattern.indexOf("[...]"));
				String postfix = filePattern.substring(filePattern.indexOf("[...]") + "[...]".length());
				
				newInput = Input.createInputSet(inputTypeMap.get(input[0]), prefix, postfix); // type and spliced name
			} else {
				newInput = Input.createInput(inputTypeMap.get(input[0]), input[1]); // type and name
			}
			inputList.add(newInput);
		}
		return inputList;
	}
	
	private Parameter parseParameter(AdvancedStringTokenizer tokens) throws ParseException {
		
		if (!"PARAMETER".equals(tokens.next())) {
			throw new ParseException("VVSADL param line should start with \"PARAMETER\", not " + tokens.current(), unitName);
		}

		String name = tokens.next();
		ParameterType type = null;
		String[] options = null;
		if (tokens.peek().startsWith("[")) {
			options = readSelectionType(tokens);
			type = ParameterType.ENUM;
		} else {
			type = ParameterType.valueOf(tokens.next());
		}
		
		String from = null;
		String to = null;
		String defaultValue = null; 

		if ("FROM".equals(tokens.peek())) {
			tokens.next(); // skip "FROM"
			from = tokens.next(); 
		}

		if ("TO".equals(tokens.peek())) {
			tokens.next(); // skip "TO"
			to = tokens.next();
		}

		if ("DEFAULT".equals(tokens.peek())) {
			tokens.next(); // skip "DEFAULT"
			defaultValue = tokens.next();
		}

		String comment = readComment(tokens);
		
		Parameter parameter = new Parameter(name, type, options, from, to, defaultValue, comment);
		
		return parameter;
	}


	private String readComment(AdvancedStringTokenizer tokens) throws ParseException {
		String comment = null;
		try {
			comment = tokens.next().substring(1); // strip (
			if (comment.contains(")")) {
				return ""; // special case:it was an empty comment
			}
			
			while (tokens.hasNext() && !tokens.peek().contains(")")) {
				comment += (" " + tokens.next());
			}
			String lastToken = tokens.next(); 
			comment += (" " + lastToken.substring(0, lastToken.length() - 1));		
			return comment;
			
		} catch (IndexOutOfBoundsException e) {
			logger.debug("on exception comment was: " + comment);
			throw new ParseException("comment start or end missing", unitName);
		}
	}
	
	private String[] readSelectionType(AdvancedStringTokenizer tokens) {
		String s = "";
		do {
			s += tokens.next();
		} while (!s.endsWith("]"));
		
		s = s.substring(1, s.length()-1); // strip "[" and "]"
		
		AdvancedStringTokenizer options = new AdvancedStringTokenizer(s, true, false, ",");
		LinkedList<String> list = new LinkedList<String>();
		while (options.hasNext()) {
			list.add(options.next());
		}
		return list.toArray(new String[0]);
	}
	
	private AdvancedStringTokenizer getTokenizer(String vvsadlString) {
		return new AdvancedStringTokenizer(vvsadlString, false, true, " \t\n\r\f/");
	}
}