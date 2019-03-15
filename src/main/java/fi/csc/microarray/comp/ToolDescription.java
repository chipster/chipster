/*
 * Created on Feb 24, 2005
 *
 */
package fi.csc.microarray.comp;

import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.description.SADLDescription.Name;
import fi.csc.microarray.description.SADLSyntax.ParameterType;


/**
 * Compute service specific versions of tool description.
 * Content is overlapping with generic SADLDescription objects, but 
 * some features are not needed here and excluded. 
 */
public class ToolDescription {

	/**
	 * Describes one parameter, such as "number of iterations".
	 * 
	 */
	public static class ParameterDescription {

		private String name;
		private String comment;
		private ParameterType type;
		
		public ParameterDescription(String name, String comment, ParameterType type) {
			this.name = name;
			this.comment = comment;
			this.type = type;
		}

		public boolean isNumeric() {
			return type.isNumeric();
		}
		
		/**
		 * Return true if this parameter is checked by the normal security policy
		 * 
		 * @return
		 */
		public boolean isChecked() {
			return this.type != ParameterType.UNCHECKED_STRING;
		}

		public String getComment() {
			return comment;
		}

		public String getName() {
			return name;
		}
	}
	
	/**
	 * Describes an output (parameter name and file name). 
	 */
	public static class OutputDescription {
        private Name fileName;
        private boolean optional;
        private boolean meta;

        public Name getFileName() {
            return fileName;
        }
        
        public boolean isOptional() {
            return optional;
        }

        public boolean isMeta() {
            return meta;
        }

	    public OutputDescription(Name fileName, boolean optional, boolean meta) {
	        this.fileName = fileName;
	        this.optional = optional;
	        this.meta = meta;
	    }
	}
	

	/**
	 * Describes an input (parameter name and file name). 
	 */
	public static class InputDescription {
        private String fileName;

        public String getFileName() {
            return fileName;
        }
	    
	    public InputDescription(String fileName) {
	        this.fileName = fileName;
	    }
	}

	
	private String id;

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
	private String displayName;
	
	/**
	 * Description.
	 */
	private String comment;

	

	private List<InputDescription> inputFiles = new LinkedList<InputDescription>();
	private List<OutputDescription> outputFiles = new LinkedList<OutputDescription>();
	private List<ParameterDescription> parameters = new LinkedList<ParameterDescription>();
	private String sourceCode;
	private String helpURL = null;

	private String initialiser;

	private Integer slotCount;
	
	public String getCommand() {
		return command;
	}
	
	public Object getImplementation() {
		return implementation;
	}

	public List<InputDescription> getInputFiles() {
		return inputFiles;
	}
	
	public List<OutputDescription> getOutputFiles() {
		return outputFiles;
	}
	
	public Iterable<ParameterDescription> getParameters() {
		return parameters;
	}
	
	public void addParameter(ParameterDescription pd) {
		parameters.add(pd);
	}
	
	public String getInitialiser() {
		return this.initialiser;
	}

	public void setInitialiser(String initialiser) {
		this.initialiser = initialiser;
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

	public String getID() {
		 return this.id;
	}
	
	public String getDisplayName() {
		if (displayName == null) {
			return id;
		} else {
			return displayName;
		}
	}

	public void setDisplayName(String displayName) {
		this.displayName = displayName;
	}
	

	public void addInputFile(String fileName) {
		inputFiles.add(new InputDescription(fileName));
	}
	
	public void addOutputFile(Name fileName, boolean optional, boolean meta) {
		outputFiles.add(new OutputDescription(fileName, optional, meta));
	}

	public void setSourceCode(String sourceCode) {
		this.sourceCode = sourceCode;		
	}

	public String getSourceCode() {
		return sourceCode;
	}

	public void setHelpURL(String helpURL) {
	    this.helpURL = helpURL;
	}

    public String getHelpURL() {
	    return helpURL;
	}
	
	public void setID(String id) {
		this.id = id;
	}

	public Integer getSlotCount() {
		if (this.slotCount == null) {
			return 1;
		} else {
			return this.slotCount;
		}
	}

	public void setSlotCount(Integer slotCount) {
		this.slotCount = slotCount;
	}
}
 