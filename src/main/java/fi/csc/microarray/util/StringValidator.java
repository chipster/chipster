package fi.csc.microarray.util;

import java.util.Locale;

import org.antlr.runtime.ANTLRStringStream;
import org.antlr.runtime.CharStream;
import org.antlr.runtime.CommonTokenStream;
import org.antlr.runtime.RecognitionException;
import org.antlr.runtime.tree.ParseTree;
import org.antlr.tool.ANTLRErrorListener;
import org.antlr.tool.ErrorManager;
import org.antlr.tool.Grammar;
import org.antlr.tool.Interpreter;
import org.antlr.tool.Message;
import org.antlr.tool.ToolMessage;

import antlr.ANTLRException;

public class StringValidator {

	private Grammar pg;
	private Grammar g;

	static {
		ErrorManager.setLocale(Locale.ENGLISH); // dont try to use localised error messages
		ErrorManager.setErrorListener(new ANTLRErrorListener() {

			public void info(String arg0) {
				System.out.println("info " + arg0);
			}

			public void error(Message arg0) {
				System.out.println("error " + arg0);
			}

			public void warning(Message arg0) {
				System.out.println("warning " + arg0);
			}

			public void error(ToolMessage arg0) {
				System.out.println("error " + arg0);
			}
			
		});
	}
	
	public StringValidator(String tokeniserGrammar, String parserGrammar) throws ANTLRException, RecognitionException {
		pg = new Grammar(parserGrammar);
		g = new Grammar();
		g.importTokenVocabulary(pg);
		g.setFileName("<string>"); // must be set or we get NPE
		g.setGrammarContent(tokeniserGrammar);
	}

	public boolean validate(String string, String startingRule) throws RecognitionException {
		CharStream input = new ANTLRStringStream(string);
		Interpreter lexEngine = new Interpreter(g, input);
		ParseTree tree;
		CommonTokenStream tokens = new CommonTokenStream(lexEngine);
		Interpreter parseEngine = new Interpreter(pg, tokens);
		tree = parseEngine.parse(startingRule);
		System.out.println(tree.toStringTree());
		return !tree.toStringTree().contains("Exception");
	}

}
