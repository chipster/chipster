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
 * <p>
 * Parses SADL analysis descriptions. SADL stands for Simple Analysis
 * Description Language. It is used to manage analysis description data inside
 * the Chipster system. Parsing is event based (cf. SAX).
 * </p>
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
		this.inputTypeMap.put(type.getName(), type);
	}

	/**
	 * Parsing is implemented with recursive descent algorithm (with 1 token
	 * look-a-head).
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

		// read possible parameters
		while (nextTokenIs(tokens, SADLSyntax.KEYWORD_PARAMETER)) {
			skip(tokens, SADLSyntax.KEYWORD_PARAMETER);
			Parameter parameter = parseParameter(tokens);
			description.addParameter(parameter);
		}

		// read possible runtime
		while (nextTokenIs(tokens, SADLSyntax.KEYWORD_RUNTIME)) {
			skip(tokens, SADLSyntax.KEYWORD_RUNTIME);
			String runtime = parseRuntime(tokens);
			description.setRuntime(runtime);
		}

		// read possible slots
		while (nextTokenIs(tokens, SADLSyntax.KEYWORD_SLOTS)) {
			skip(tokens, SADLSyntax.KEYWORD_SLOTS);
			int slotCount = parseInt(tokens);
			description.setSlotCount(slotCount);
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
			name.setPostfix(rawName.substring(
					rawName.indexOf(SADLSyntax.NAME_SET_DESIGNATOR) + SADLSyntax.NAME_SET_DESIGNATOR.length()));

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

	private String parseRuntime(SADLTokeniser tokens) throws ParseException {
		return tokens.next();
	}

	private int parseInt(SADLTokeniser tokens) throws ParseException {
		String nextToken = tokens.next();
		try {
			return Integer.parseInt(nextToken);
		} catch (NumberFormatException e) {
			throw new ParseException("slot count is not integer: " + nextToken);
		}
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
			this.validateFromOrToOrDefaultAsNumber(type, from, name, SADLSyntax.KEYWORD_FROM);
		}

		if (nextTokenIs(tokens, SADLSyntax.KEYWORD_TO)) {
			skip(tokens, SADLSyntax.KEYWORD_TO);
			to = tokens.next();
			this.validateFromOrToOrDefaultAsNumber(type, to, name, SADLSyntax.KEYWORD_TO);
		}

		// validate from and to relative to each other
		// from and to have been validated as numbers before
		if (from != null && to != null) {
			double fromDouble = Double.parseDouble(from);
			double toDouble = Double.parseDouble(to);
			if (toDouble < fromDouble) {
				throw new ParseException("parameter " + name.getID() + " " + SADLSyntax.KEYWORD_TO + " '" + to
						+ "' is less than " + SADLSyntax.KEYWORD_FROM + " '" + from + "'");
			}
		}

		if (nextTokenIs(tokens, SADLSyntax.KEYWORD_DEFAULT)) {
			skip(tokens, SADLSyntax.KEYWORD_DEFAULT);
			defaultValues = parseDefaultValues(tokens);

			// validate default against from and to
			if (type == ParameterType.DECIMAL || type == ParameterType.INTEGER) {
				this.validateDefaultNumber(defaultValues, name, type, from, to);
			}
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

	private void validateFromOrToOrDefaultAsNumber(ParameterType type, String value, Name name,
			String fromOrToOrDefault) throws ParseException {
		if (type == ParameterType.INTEGER) {
			try {
				Integer.parseInt(value);
			} catch (NumberFormatException e) {
				throw new ParseException("parameter " + name.getID() + " has illegal " + fromOrToOrDefault + " value '"
						+ value + "', should be an integer");
			}
		} else if (type == ParameterType.DECIMAL) {
			try {
				Double.parseDouble(value);
			} catch (NumberFormatException e) {
				throw new ParseException("parameter " + name.getID() + " has illegal " + fromOrToOrDefault + " value '"
						+ value + "', should be a decimal number");
			}
		}
	}

	private void validateDefaultNumber(String[] defaultValues, Name name, ParameterType type, String from, String to)
			throws ParseException {
		// only one default for INTEGER and DECIMAL
		if (defaultValues.length != 1) {
			throw new ParseException(
					"parameter " + name.getID() + " has illegal number of defaults '" + defaultValues.length
							+ "', should be '1' for " + ParameterType.INTEGER + " " + ParameterType.DECIMAL);
		}
		String defaultString = defaultValues[0];

		// parses as an INTEGER or DECIMAL
		this.validateFromOrToOrDefaultAsNumber(type, defaultString, name, SADLSyntax.KEYWORD_DEFAULT);

		// not less than FROM
		if (from != null && Double.parseDouble(from) > Double.parseDouble(defaultString)) {
			throw new ParseException("parameter " + name.getID() + " " + SADLSyntax.KEYWORD_DEFAULT + " '"
					+ defaultString + "' is less than " + SADLSyntax.KEYWORD_FROM + " '" + from + "'");
		}

		// not greater than TO
		if (to != null && Double.parseDouble(to) < Double.parseDouble(defaultString)) {
			throw new ParseException("parameter " + name.getID() + " " + SADLSyntax.KEYWORD_DEFAULT + " '"
					+ defaultString + "' is greater than " + SADLSyntax.KEYWORD_TO + " '" + to + "'");
		}

	}

}