package fi.csc.microarray.description;

import fi.csc.microarray.databeans.Dataset;


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
 * Quoted texts are snippets of SADL. For example, TOOL is a term
 * that is rewritten using the given rules, but "TOOL" is a string that
 * should be found in the source code. Operators ?, +, * and | have their
 * common semantics.</p>
 * 
 * <pre>
 * -> TOOL+
 * 
 * TOOL -> "TOOL" CATEGORY_NAME "/" NAME DESCRIPTION INPUT* METAINPUT* OUTPUT* METAOUTPUT* PARAMETER*
 * CATEGORY_NAME -> TOKEN
 * INPUT -> "INPUT" INPUT_DESCRIPTION
 * METAINPUT -> "METAINPUT" INPUT_DESCRIPTION
 * INPUT_DESCRIPTION -> OPTIONALITY? NAME "TYPE" TYPE_NAME
 * TYPE_NAME -> TOKEN (see SADLSyntax.ParameterType for declaration, implementations pluggable)
 * OUTPUT -> "OUTPUT" OUTPUT_DESCRIPTION
 * METAOUTPUT -> "METAOUTPUT" OUTPUT_DESCRIPTION
 * OUTPUT_DESCRIPTION -> OPTIONALITY? NAME 
 * PARAMETER -> "PARAMETER" OPTIONALITY? NAME "TYPE" PARAMETER_TYPE PARAMETER_FROM? PARAMETER_TO? PARAMETER_DEFAULT? DESCRIPTION 
 * PARAMETER_TYPE -> TOKEN | PARAMETER_TYPE_ENUM
 * PARAMETER_TYPE_ENUM -> "[" PARAMETER_TYPE_ENUM_ELEMENTS "]"
 * PARAMETER_TYPE_ENUM_ELEMENTS -> NAME | NAME "," PARAMETER_TYPE_ENUM_ELEMENTS
 * PARAMETER_FROM -> "FROM" TOKEN
 * PARAMETER_TO -> "TO" TOKEN
 * PARAMETER_DEFAULT -> "DEFAULT" PARAMETER_DEFAULT_ELEMENT
 * PARAMETER_DEFAULT_ELEMENT -> TOKEN | TOKEN "," PARAMETER_DEFAULT_ELEMENT 
 * OPTIONALITY -> "OPTIONAL"
 * NAME -> TOKEN | TOKEN ":" TOKEN
 * DESCRIPTION -> TOKEN
 * TOKEN -> any single token produced by tokeniser
 * </pre>
 * 
 * <p>There is special handling for the NAME tokens: they are tokenized normally, but parser emits a different type of name 
 * if the token contains SADLSyntax.NAME_SET_DESIGNATOR (meaning that we are dealing with name set instead of single name).</p>
 * 
 * <p>TOKEN refers to any single token produced by tokeniser (SADLTokeniser). Tokens can be keywords, operators, 
 * strings, quoted string or strings in parentheses. String in parentheses are strongly recommended for descriptions.</p>  
 *  
 * <p>Below is an example of a SADL snippet.</p>
 * <pre>
 * TOOL "Test utilities" / util-test.R: "Test tool" (Just a test analysis for development. These descriptions are sometimes very
 * long and might get hard to read.)
 * INPUT microarray{...}.tsv: "Raw data files" TYPE CDNA
 * INPUT phenodata.tsv: "Experiment description" TYPE GENERIC
 * OUTPUT result{...}.txt: "Result files"
 * OUTPUT OPTIONAL error.txt: "Error, if any"
 * PARAMETER value1: "The first value" TYPE INTEGER FROM 0 TO 200 DEFAULT 10 (the first value of the result set)
 * PARAMETER value2: "The second value" TYPE DECIMAL FROM 0 TO 200 DEFAULT 20 (the second value of the result set)
 * PARAMETER OPTIONAL value3: "The third value" TYPE DECIMAL FROM 0 TO 200 DEFAULT 30.2 (the third value of the result set)
 * PARAMETER value4: "The fourth value" TYPE PERCENT DEFAULT 34 (how much we need)
 * PARAMETER method: "The enumeration" TYPE [option1: "First option", option2: "Second option", option3: "Third option"] FROM 1 TO 2 DEFAULT option1, option2 (which options are selected)
 * PARAMETER genename: "Gene name" TYPE STRING DEFAULT at_something (which gene we are interested in)
 * PARAMETER key: "Key column" TYPE COLUMN_SEL (which column we use as a key)
 * </pre>  
 * 
 * @see fi.csc.microarray.description.SADLTokeniser
 * @see fi.csc.microarray.description.SADLParser
 * 
 * @author Aleksi Kallio
 *
 */
public class SADLSyntax {

	public static final String KEYWORD_DEFAULT = "DEFAULT";
	public static final String KEYWORD_TO = "TO";
	public static final String KEYWORD_FROM = "FROM";
	public static final String KEYWORD_OPTIONAL = "OPTIONAL";
	public static final String KEYWORD_TYPE = "TYPE";
	public static final String KEYWORD_PARAMETER = "PARAMETER";
	public static final String KEYWORD_METAOUTPUT = "METAOUTPUT";
	public static final String KEYWORD_OUTPUT = "OUTPUT";
	public static final String KEYWORD_METAINPUT = "METAINPUT";
	public static final String KEYWORD_INPUT = "INPUT";
	public static final String KEYWORD_TOOL = "TOOL";

	public static final String NAME_SET_DESIGNATOR = "{...}";
	public static final String NAME_SEPARATOR = ":";
	public static final String CATEGORY_SEPARATOR = "/";
	public static final String ENUM_OPEN = "[";
	public static final String ENUM_CLOSE = "]";
	public static final String LIST_SEPARATOR = ",";
	public static final String COMMENT_OPEN = "(";
	public static final String COMMENT_CLOSE = ")";
	public static final String QUOTE = "\"";
	public static final String ESCAPE = "\\";
	
	public static interface InputType {
		public boolean isTypeOf(Dataset dataBean);
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
