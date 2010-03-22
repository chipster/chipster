package fi.csc.microarray.description;

import fi.csc.microarray.databeans.DataBean;


/**
 * <p>SADL (Simple Analysis Description Language) is a simple language for 
 * describing analysis operations so that they can be used in the Chipster system. Operations
 * have names and are categorised. They have both inputs and parameters. Inputs are data
 * files (large) and parameters are values (small) that user gives when calling the 
 * operation. Operations produce outputs (data files) and output text i.e. what they write
 * to sysout.</p>
 *    
 * <p>Syntax defination is below. It is in the form of rewrite rules.
 * First rule in the list is the initial rule where rewriting is started. 
 * Hyphenated texts are snippets of SADL. For example, TOOL is a term
 * that is rewritten using the given rules, but "TOOL" is a string that
 * should be found in the source code. Operators ?, +, * and | have their
 * common semantics. Lower case identifiers: number, string and string_no_ws (no 
 * whitespace) should be obvious.</p>
 * 
 * <pre>
 * -> TOOL+
 * TOOL -> "TOOL" NAME "/" NAME DESCRIPTION INPUT? OUTPUT? PARAMETER*
 * NAME -> NAME_SINGLE | NAME_SET  
 * NAME_SINGLE -> string_no_ws | """ string """
 * NAME_SET -> string_no_ws "[...]" string_no_ws
 * INPUT -> "INPUT" INPUT_TYPE NAME | INPUT_TYPE FILENAMESET
 * INPUTTYPE -> [see SADLSyntax.InputType for declaration, implementations pluggable]
 * OUTPUT -> "OUTPUT" FILENAME
 * PARAMETER -> "PARAMETER" NAME PARAMTYPE RANGE? DEFAULT? DESCRIPTION 
 * PARAMTYPE -> [see SADLSyntax.ParameterType for declaration, implementations pluggable]
 * RANGE -> "FROM" number "TO" number
 * DEFAULT -> "DEFAULT" number
 * DESCRIPTION -> "(" string ") 
 * </pre>  
 * 
 * <p>Below is an example of a SADL snippet.</p>
 * <pre>
 * ANALYSIS Test/test (Just a test analysis for development)
 * INPUT CDNA microarray[...].txt OUTPUT results.txt, messages.txt
 * PARAMETER value1 INTEGER FROM 0 TO 200 DEFAULT 10 (the first value of the result set)
 * PARAMETER value2 DECIMAL FROM 0 TO 200 DEFAULT 20 (the second value of the result set)
 * PARAMETER value3 DECIMAL FROM 0 TO 200 DEFAULT 30.2 (the third value of the result set)
 * PARAMETER method PERCENT DEFAULT 34 (how much we need)
 * PARAMETER method [linear, logarithmic, exponential] DEFAULT logarithmic (which method to apply)
 * PARAMETER genename STRING DEFAULT at_something (which gene we are interested in)
 * PARAMETER key COLNAME (which column we use as a key)
 * </pre>  
 * 
 * @see fi.csc.microarray.description.SADLParser
 * @see SADLSyntax.InputType 
 * @see SADLSyntax.ParameterType
 * 
 * @author Aleksi Kallio
 *
 */
public class SADLSyntax {
	
	public static final String INPUT_SET_DESIGNATOR = "{...}";
	public static final String NAME_SEPARATOR = ":";
	public static final String CATEGORY_SEPARATOR = "/";
	
	public static interface InputType {
		public boolean isTypeOf(DataBean dataBean);
		public String getName();
		public boolean isMetadata();
	}
	
	public static enum ParameterType {
		/**
		 * Integer number.
		 */
		INTEGER,
		
		/**
		 * Decimal number.
		 */
		DECIMAL,
		
		/**
		 * Integer between 0 and 100 (inclusive).
		 */
		PERCENT,
		
		/**
		 * A character string.
		 */
		STRING,
		
		/**
		 * Enumeration from a set of given values (specified as this type is referred). 
		 */
		ENUM,
		
		/**
		 * Name of input matrix column, for choosing columns from inputs.
		 */
		COLUMN_SEL,
		
		/**
		 * Name of metainput matrix column, for choosing columns from metainputs.
		 */
		METACOLUMN_SEL,
		
		/**
		 * Name of input, for choosing from multiple input datasets.
		 */
		INPUT_SEL;

		public static boolean isValid(String typeName) {
			for (ParameterType type : values()) {
				if (type.toString().equals(typeName)) {
					return true;
				}
			}
			return false;
		}
		
		public boolean isNumeric() {
			return this == INTEGER || this == DECIMAL || this == PERCENT;
		}

	}
}
