package fi.csc.microarray.cluster;

import fi.csc.microarray.util.LookaheadStringReader;

public class TreeParseException extends Exception {

	public TreeParseException(String desc, LookaheadStringReader reader) {
		super(desc + " (around ..." + reader.around(10) + "...)");
	}
}
