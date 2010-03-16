package fi.csc.microarray.description;

import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.description.SADLParser.ParseException;

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

	private static class CharTokeniser {

		private String sadl;
		private int index = 0;
		private String unitName;

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
				if (peek() == '(') {

					next(); // skip '('
					
					while (!atEnd() && peek() != ')') {
						token += next();
					}

					next(); // skip ')'

				// read special quoted block
				} else if (peek() == '"') {

					next(); // skip '"'

					while (!atEnd() && peek() != '"') {
						token += next();
					}

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

		private boolean atEnd() {
			return index >= sadl.length();
		}

		private boolean isWhiteSpace() {
			char c = peek();
			return Character.isWhitespace(c);
		}

		private boolean isOperator() {
			char c = peek();
			return c == ':' || c == '[' || c == ']' || c == ',' || c == '/';
		}

		public char next() {
			return sadl.charAt(index++);
		}

		public char peek() {
			return sadl.charAt(index);
		}
	}
}
