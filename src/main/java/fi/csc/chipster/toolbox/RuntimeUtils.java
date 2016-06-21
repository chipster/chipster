package fi.csc.chipster.toolbox;

public class RuntimeUtils {
	/**
	 * Use runtime name to figure out the dir where tools
	 * script file is located.
	 *  
	 */
	public static String getToolDirFromRuntimeName(String runtime) {
		if (runtime == null || runtime.isEmpty()) {
			throw new IllegalArgumentException(runtime);
		} else if (runtime.startsWith("python")) {
			return "python";
		} else if (runtime.startsWith("java")) {
			return "java";

			// add non-R stuff starting with R before this
		} else if (runtime.startsWith("R")) {
			return "R";
		} else {
			return runtime;
		}
	}

}
