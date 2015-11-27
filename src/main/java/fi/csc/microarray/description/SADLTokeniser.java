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

	public static enum TokenType {
		DESCRIPTION,
		QUOTED,
		OPERATOR,
		NORMAL
	}
	
	private List<String> tokens;
	private List<TokenType> types;
	int index = -1;

	public SADLTokeniser(String sadl, String unitName) throws ParseException {
		CharTokeniser t = new CharTokeniser(sadl, unitName);
		t.tokenise();
		this.tokens = t.getTokens();
		this.types = t.getTypes();
		
		if (this.tokens.size() != this.types.size()) {
			throw new RuntimeException("internal parsing error, number of types and tokens does not match");
		}
	}
	
	public String peek() {
		if (!hasNext()) {
			return null;
		}
		return tokens.get(index + 1);
	}

	public TokenType peekType() {
		if (!hasNext()) {
			return null;
		}
		return types.get(index + 1);
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
		private int charIndex = 0;
		private String unitName;
		private boolean escaped;
		private boolean escapingEnabled;
		private List<String> tokens = new LinkedList<String>();
		private List<TokenType> types = new LinkedList<TokenType>();
		
		public CharTokeniser(String sadl, String unitName) {
			this.sadl = sadl;
			this.unitName = unitName;
		}

		public void tokenise() throws ParseException {

			// iterate through the whole string
			while (!atEnd()) {

				int startedAt = charIndex;
				String token = "";
				TokenType type = TokenType.NORMAL; 

				// read special comment block
				if (peek() == SADLSyntax.COMMENT_OPEN.charAt(0)) {

					next(); // skip '('
					setEscapingEnabled(true);
					
					while (!atEnd() && !isOperator(SADLSyntax.COMMENT_CLOSE.charAt(0))) {
						token += next();
					}

					setEscapingEnabled(false);
					next(); // skip ')'
					
					type = TokenType.DESCRIPTION;

				// read special quoted block
				} else if (peek() == SADLSyntax.QUOTE.charAt(0)) {

					next(); // skip '"'
					setEscapingEnabled(true);

					while (!atEnd() && !isOperator(SADLSyntax.QUOTE.charAt(0))) {
						token += next();
					}

					setEscapingEnabled(false);
					next(); // skip '"'

					type = TokenType.QUOTED;

				// read special operator block
				} else if (isOperator()) {

					while (!atEnd() && isOperator()) {
						token += next();
					}
					
					type = TokenType.OPERATOR;

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
				if (charIndex == startedAt) {
					throw new ParseException("tokeniser got stuck at ", unitName);
				}

				// emit token
				tokens.add(token);
				types.add(type);

			}
		}

		public List<String> getTokens() {
			return tokens;
		}

		public List<TokenType> getTypes() {
			return types;
		}

		private void setEscapingEnabled(boolean b) {
			this.escapingEnabled = true;
		}

		private boolean atEnd() {
			return charIndex >= sadl.length();
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
			charIndex++;
			return c;
		}

		public char peek() {
			
			// skip escape char and remember the we have escaped the next char
			this.escaped = false;
			if (escapingEnabled && sadl.charAt(charIndex) == SADLSyntax.ESCAPE.charAt(0)) {
				charIndex++;
				escaped = true;
			}
			
			return sadl.charAt(charIndex);
		}
	}
}
