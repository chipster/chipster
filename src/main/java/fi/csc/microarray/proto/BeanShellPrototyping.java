package fi.csc.microarray.proto;

import java.io.InputStreamReader;

import bsh.EvalError;
import bsh.Interpreter;

public class BeanShellPrototyping {
	
	public static void main(String[] args) throws EvalError  {
		Interpreter i = new Interpreter();
		i.set("h", 5);
		i.eval(new InputStreamReader(BeanShellPrototyping.class.getResourceAsStream("/bsh/test.bsh")));		
	}

}
