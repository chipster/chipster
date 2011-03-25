package fi.csc.microarray.description;

import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.description.SADLParser.ParseException;

/**
 * Parses string into tokens, suitable for SADLParser. The tokeniser has some understating
 * of SADL syntax so that it can produce higher level tokens and make parsing easier.
 * Tokens can be keywords, operators, strings, quoted string or strings in parentheses.
 *
 * @see SADLParser
 * 
 * @author Aleksi Kallio
 *
 */
public class SADLTokeniser {

	private List<String> tokens;
	int index = -1;

	public SADLTokeniser(String sadl, String unitName) throws ParseException {
		CharTokeniser t = new CharTokeniser(sadl, unitName);
		this.tokens = t.tokenise();
	}
	
	public String peek() {
		return tokens.get(index + 1);
	}

	public boolean hasNext() {
		return (index+1) < tokens.size();
	}

	public String next() {
		if (!hasNext()) {
			throw new RuntimeException("tokeniser error: no tokens left (there were " + tokens.size() + " tokens)");
		}
		return tokens.get(++index);
	}

	public static String[] blockEndingOperators() {
		return new String[] {
				SADLSyntax.QUOTE,
				SADLSyntax.COMMENT_CLOSE
		};
	}
	
	public static String[] tokenEndingOperators() {
		return new String[] {
				SADLSyntax.NAME_SEPARATOR,
				SADLSyntax.ENUM_OPEN,
				SADLSyntax.ENUM_CLOSE,
				SADLSyntax.LIST_SEPARATOR
		};
	}

	private static class CharTokeniser {

		private String sadl;
		private int index = 0;
		private String unitName;
		private boolean escaped;
		private boolean escapingEnabled;

		public CharTokeniser(String sadl, String unitName) {
			this.sadl = sadl;
			this.unitName = unitName;
		}

		public List<String> tokenise() throws ParseException {

			List<String> tokens = new LinkedList<String>();

			// iterate through the whole string
			while (!atEnd()) {

				int startedAt = index;
				String token = "";

				// read special comment block
				if (peek() == SADLSyntax.COMMENT_OPEN.charAt(0)) {

					next(); // skip '('
					setEscapingEnabled(true);
					
					while (!atEnd() && !isOperator(SADLSyntax.COMMENT_CLOSE.charAt(0))) {
						token += next();
					}

					setEscapingEnabled(false);
					next(); // skip ')'

				// read special quoted block
				} else if (peek() == SADLSyntax.QUOTE.charAt(0)) {

					next(); // skip '"'
					setEscapingEnabled(true);

					while (!atEnd() && !isOperator(SADLSyntax.QUOTE.charAt(0))) {
						token += next();
					}

					setEscapingEnabled(false);
					next(); // skip '"'

				// read special operator block
				} else if (isOperator()) {

					while (!atEnd() && isOperator()) {
						token += next();
					}

				// read regular block of non-whitespace
				} else {

					while (!atEnd() && !isWhiteSpace() && !isOperator()) {
						token += next();
					}
				}

				// read following block of whitespace
				while (!atEnd() && isWhiteSpace()) {
					next();
				}

				// check that we made progress
				if (index == startedAt) {
					throw new ParseException("tokeniser got stuck at ", unitName);
				}

				// emit token
				tokens.add(token);

			}

			return tokens;
		}

		private void setEscapingEnabled(boolean b) {
			this.escapingEnabled = true;
		}

		private boolean atEnd() {
			return index >= sadl.length();
		}

		private boolean isWhiteSpace() {
			char c = peek();
			return Character.isWhitespace(c);
		}

		private boolean isOperator(char operator) {
			boolean isOperator =  peek() == operator;
			if (escaped) {
				return false; // escaped chars are never interpreted as operators
			}
			return isOperator;
		}
		
		private boolean isOperator() {
			char c = peek();
			if (escaped) {
				return false; // escaped chars are never interpreted as operators
			}
			
			for (String operator : tokenEndingOperators()) {
				if (c == operator.charAt(0)) {
					return true;
				}
			}
			return false;
		}

		public char next() {
			char c = peek();
			index++;
			return c;
		}

		public char peek() {
			
			// skip escape char and remember the we have escaped the next char
			this.escaped = false;
			if (escapingEnabled && sadl.charAt(index) == SADLSyntax.ESCAPE.charAt(0)) {
				index++;
				escaped = true;
			}
			
			return sadl.charAt(index);
		}
	}
}
