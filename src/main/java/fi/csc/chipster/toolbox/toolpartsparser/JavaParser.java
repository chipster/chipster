package fi.csc.chipster.toolbox.toolpartsparser;

import java.io.File;

import org.apache.log4j.Logger;

import fi.csc.chipster.toolbox.SADLTool.ParsedScript;
import fi.csc.microarray.comp.java.JavaCompJobBase;

public class JavaParser implements ToolPartsParser {

	static final Logger logger = Logger.getLogger(JavaParser.class);


	@Override
	public ParsedScript parse(File moduleDir, String sourceResourceName) throws ClassNotFoundException, InstantiationException, IllegalAccessException {
		// get the job class
		Class<? extends Object> jobClass = null;

		jobClass = Class.forName(sourceResourceName);
		
		assert(JavaCompJobBase.class.isAssignableFrom(jobClass));
		JavaCompJobBase jobInstance; 
		jobInstance = (JavaCompJobBase)jobClass.newInstance();
		
		// TODO what to do with other parts
		ParsedScript ps = new ParsedScript();
		ps.SADL = jobInstance.getSADL();
		ps.source = "Source code for this tool is available within Chipster source code.";
		
		return ps;
	}
}
