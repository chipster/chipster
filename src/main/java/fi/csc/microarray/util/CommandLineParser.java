package fi.csc.microarray.util;

import java.util.LinkedHashMap;
import java.util.Map;

public class CommandLineParser {

	public static final String NON_VALUE = "";
	
	public static class CommandLineException extends Exception {
		public CommandLineException(String message) {
			super(message);
		}
	}
	
	private static class Parameter {
		String name;
		String value = null;		
		String defaultValue = null;
		String description = null;
		boolean mandatory = false;
		boolean needsValue = false;
	}
	
	private Map<String, CommandLineParser.Parameter> parameters = new LinkedHashMap<String, CommandLineParser.Parameter>();
	private boolean userAskedHelp = false;
	
	public void addParameter(String name, boolean mandatory, boolean hasValue, String defaultValue, String description) {
		Parameter p = new Parameter();
		p.name = name;
		p.value = null;
		p.needsValue = hasValue;
		p.mandatory = mandatory;
		p.defaultValue = defaultValue;
		p.description = description;
		
		parameters.put(p.name, p);
	}
	
	public void parse(String[] args) throws CommandLineException {
		for (int i = 0; i < args.length; i++) {
			
			if ("help".equals(args[i]) || 
					"-help".equals(args[i]) || 
					"--help".equals(args[i]) || 
					"-h".equals(args[i])) {
				userAskedHelp = true;
				return;
			}
			
			Parameter parameter = parameters.get(args[i]);
			if (parameter == null) {
				throw new CommandLineException("unknown parameter: " +  args[i]);
			}
			if (parameter.needsValue) {
				try {
					parameter.value = args[i+1];
					i++; // do not treat value as a new parameter
				} catch (ArrayIndexOutOfBoundsException e) {
					throw new CommandLineException("parameter " +  parameter.name + " is missing a value");
				}
			} else {
				parameter.value = NON_VALUE; 
			}
		}
		
		for (Parameter p : parameters.values()) {
			if (p.mandatory && p.value == null) {
				throw new CommandLineException("parameter " + p.name + " is mandatory, but missing");
			}
		}

		// add -help
		Parameter help = new Parameter();
		help.name = "-help";
		help.description = "get help on command line syntax (this text)";
		parameters.put(help.name, help);	
	}
	
	public boolean hasValue(String name) throws CommandLineException {
		return getValue(name) != null;
	}
	
	public String getValue(String name) throws CommandLineException {
		Parameter p = parameters.get(name);
		if (p == null) {
			throw new CommandLineException("no such parameter: " + name);
		}
		if (p.value == null && p.defaultValue != null) {
			return p.defaultValue;
		}		
		return p.value;
	}
	
	public boolean userAskedHelp() {
		return userAskedHelp;
	}
	
	public String getDescription() {
		String desc = "";
		for (Parameter p : parameters.values()) {
			String line = "";
			String header = p.needsValue ? (p.name + " <value>") : p.name; 
			int padding = 8 - header.length();
			String pad = Strings.repeat(" ", padding > 0 ? padding : 0); 
			line += header + pad;
			line += p.mandatory ? "mandatory" : "";
			line += (p.defaultValue != null) ? (", default:" + p.defaultValue) : "";
			int tabs = (int)Math.ceil(((4*8) - line.length()) / 8.0f);
			line += (p.description != null) ? (Strings.repeat("\t", tabs) + p.description) : "";
			line += "\n";
			desc += line;
		}
		return desc;
	}

}
