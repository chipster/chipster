/*
 * Created on Feb 24, 2005
 *
 */
package fi.csc.microarray.analyser;

import java.util.Date;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.message.JobMessage;


/**
 * Describes one analysis, such as "median normalisation".
 * 
 * @author hupponen, akallio 
 */
public class AnalysisDescription {
	/**
	 * Logger for this class
	 */
	@SuppressWarnings("unused")
	private static final Logger logger = Logger
			.getLogger(AnalysisDescription.class);

	/**
	 * Describes one parameter, such as "number of iterations".
	 * 
	 * @author akallio
	 */
	public static class ParameterDescription {

		private String name;
		private String comment;
		private boolean numeric;
		
		public ParameterDescription(String name, String comment, boolean numeric) {
			this.name = name;
			this.comment = comment;
			this.numeric = numeric;
		}

		public boolean isNumeric() {
			return numeric;
		}

		public String getComment() {
			return comment;
		}

		public String getName() {
			return name;
		}
	}
	
	/**
	 * Actual executable that handles the analysis.
	 */
	private String command;
	
	/**
	 * The actual content of the operation implementation.
	 */
	private Object implementation;
	
	/**
	 * Analysis name (used in GUI etc.)
	 */
	private String name;
	
	/**
	 * Description.
	 */
	private String comment;

	
	
	private List<String> outputFiles = new LinkedList<String>();
	private List<ParameterDescription> parameters = new LinkedList<ParameterDescription>();
	private String sourceCode;
	private String category;
	private String vvsadl;
	private AnalysisHandler handler;
	
	// these are needed for update check
	/** Name of the original source script or java class or whatever */
	private String sourceResourceName;
	private String sourceResourceFullPath;
	
	private Date creationTime = new Date();

	/**
	 * Initializes empty (non-usable) description.
	 *
	 */
	public AnalysisDescription(AnalysisHandler handler) {
		this.handler = handler;
	}

	public String getCommand() {
		return command;
	}
	
	/**
	 * Factory method, infers correct job type.
	 * @param message
	 * @param resultHandler
	 * @return
	 * @throws AnalysisException 
	 */
	public AnalysisJob createAnalysisJob(JobMessage message, ResultCallback resultHandler) throws AnalysisException {
		return handler.createAnalysisJob(message, this, resultHandler);		
	}
	
	public Object getImplementation() {
		return implementation;
	}
	
	public List<String> getOutputFiles() {
		return outputFiles;
	}
	
	public Iterable<ParameterDescription> getParameters() {
		return parameters;
	}
	
	public void addParameter(ParameterDescription pd) {
		parameters.add(pd);
	}
	
	public static final String getStaticInitialiser() {
		return "";
	}


	public void setCommand(String command) {
		this.command = command;
	}

	public void setImplementation(Object implementation) {
		this.implementation = implementation;
	}

	public String getComment() {
		return comment;
	}

	public void setComment(String comment) {
		this.comment = comment;
	}

	public String getFullName() {
		return "\"" + getCategory() + "\"/\"" + getName() + "\""; 
	}
	
	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public void addOutputFile(String file) {
		outputFiles.add(file);
	}

	public void setSourceCode(String sourceCode) {
		this.sourceCode = sourceCode;		
	}

	public String getSourceCode() {
		return sourceCode;
	}

	public void setCategory(String category) {
		this.category = category;		
	}

	public String getCategory() {
		return category;
	}

	public void setVVSADL(String vvsadl) {
		this.vvsadl = vvsadl;
	}
	
	public String getVVSADL() {
		
		if (vvsadl != null) {
			return vvsadl;
		} else {
			return "ANALYSIS " + getFullName() + " (" + getComment() + ")";
		}
	}

	public String getSourceResourceName() {
		return sourceResourceName;
	}

	public void setSourceResourceName(String sourceResourceName) {
		this.sourceResourceName = sourceResourceName;
	}

	public String getSourceResourceFullPath() {
		return sourceResourceFullPath;
	}

	public void setSourceResourceFullPath(String sourceResourceFullPath) {
		this.sourceResourceFullPath = sourceResourceFullPath;
	}

	
	public long getCreationTime() {
		return this.creationTime.getTime();
	}
	
	public boolean isUptodate() {
		return this.handler.isUptodate(this);
	}
	
	public AnalysisHandler getHandler() {
		return this.handler;
	}

}
 