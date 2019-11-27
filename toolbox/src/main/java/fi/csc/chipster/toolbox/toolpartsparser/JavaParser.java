package fi.csc.chipster.toolbox.toolpartsparser;

import java.lang.reflect.InvocationTargetException;
import java.nio.file.Path;

import org.apache.log4j.Logger;

import fi.csc.chipster.toolbox.SADLTool.ParsedScript;
import fi.csc.microarray.comp.java.JavaCompJobBase;

public class JavaParser implements ToolPartsParser {

	static final Logger logger = Logger.getLogger(JavaParser.class);


	@Override
	public ParsedScript parse(Path moduleDir, String toolId) throws ClassNotFoundException, InstantiationException, IllegalAccessException, IllegalArgumentException, InvocationTargetException, NoSuchMethodException, SecurityException {

		String sourceResourceName = toolId.substring(0, toolId.lastIndexOf(".java")); 

		// get the job class
		Class<? extends Object> jobClass = null;

		jobClass = Class.forName(sourceResourceName);
		
		assert(JavaCompJobBase.class.isAssignableFrom(jobClass));
		JavaCompJobBase jobInstance; 
		jobInstance = (JavaCompJobBase)jobClass.getDeclaredConstructor().newInstance();
		
		// TODO what to do with other parts
		ParsedScript ps = new ParsedScript();
		ps.SADL = jobInstance.getSADL();
		ps.source = "Source code for this tool is available within Chipster source code.";
		
		return ps;
	}
}
