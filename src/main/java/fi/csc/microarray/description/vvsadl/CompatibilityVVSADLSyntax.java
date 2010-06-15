package fi.csc.microarray.description.vvsadl;

import fi.csc.microarray.databeans.Dataset;


/**
 * <p>VVSADL (Very very simple analysis description language) is a simple language for
 * describing analysis operations so that they can be used in the NAMI system. Operations
 * have names and are categorised. They have both inputs and parameters. Inputs are data
 * files (large) and parameters are values (small) that user gives when calling the
 * operation. Operations produce outputs (data files) and output text ie. what they write
 * to sysout.</p>
 *
 * <p>Syntax defination is below. It is in the form of rewrite rules.
 * First rule in the list is the initial rule where rewriting is started.
 * Quoted texts are snippets of VVSADL. For example, ANALYSIS is a term
 * that is rewritten using the given rules, but "ANALYSIS" is a string that
 * should be found in the source code. Operators ?, +, * and | have their
 * common semantics. Lower case identifiers: number, string and string_no_ws (no
 * whitespace) should be obvious.</p>
 *
 * <pre>
 * -> ANALYSIS+
 * ANALYSIS -> "ANALYSIS" CATEGORY "/" OPNAME DESCRIPTION INPUT? OUTPUT? PARAMETER*
 * CATEGORY -> NAMESTRING
 * OPNAME -> NAMESTRING
 * NAMESTRING -> string_no_ws | """ string """
 * INPUT -> INPUTTYPE FILENAME | INPUTTYPE FILENAMESET
 * INPUTTYPE -> [see VVSADLSyntax.InputType for declaration, implementations pluggable]
 * OUTPUT -> FILENAME
 * FILENAME -> string_no_ws
 * FILENAMESET -> string_no_ws "[...]" string_no_ws
 * PARAMETER -> NAME PARAMTYPE RANGE? DEFAULT? DESCRIPTION
 * PARAMTYPE -> [see VVSADLSyntax.ParameterType for declaration, implementations pluggable]
 * RANGE -> "FROM" number "TO" number
 * DEFAULT -> "DEFAULT" number
 * DESCRIPTION -> "(" string ")
 * </pre>
 *
 * <p>Below is an example of a VVSADL snippet.</p>
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
 * @see fi.csc.microarray.description.VVSADLParser
 * @see VVSADLSyntax.InputType
 * @see VVSADLSyntax.ParameterType
 *
 * @author Aleksi Kallio
 *
 */
public class CompatibilityVVSADLSyntax {


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
