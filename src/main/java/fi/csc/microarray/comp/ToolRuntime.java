package fi.csc.microarray.comp;



/**
 * Contains information about a single runtime, specified in the runtime.xml.
 * 
 * A runtime specifies an environment, in which a tool (often a script file)
 * is actually run.Examples of runtimes include different versions of R or
 *  the beanshell runtime for running BeanShell scripts.
 *  
 *  In practice, a runtime configuration consists of a name and handler
 *	configuration. Runtime name is used as a reference in the tools.xml and
 *	the handler element specifies the class which takes care of actually 
 *	running the tool.
 * 
 */
public class ToolRuntime {

	private String name;
	private JobFactory jobFactory;
	private boolean disabled;
	
	/**
	 * Setting a runtime disabled, prevents tools from being run in the specific 
	 * runtime. The disabled/enabled state of the runtime reflects the 
	 * state of the handler, which actually takes care of running the tools 
	 * 
	 * The disabled state can be specified in the runtimes.xml.
	 * 
	 * @param name
	 * @param jobFactory
	 * @param disabled
	 */
	public ToolRuntime(String name, JobFactory jobFactory, boolean disabled) {
		this.name = name;
		this.jobFactory = jobFactory;
		this.disabled = disabled;
	}

	public String getName() {
		return name;
	}

	public JobFactory getJobFactory() {
		return jobFactory;
	}

	public boolean isDisabled() {
		return disabled || jobFactory.isDisabled();
	}

}
