package fi.csc.microarray.proto.repository;

/**
 * Partial implementation that returns null when feature is not required.
 * 
 * @author akallio
 */
public abstract class ExperimentBase implements Experiment {

	public String getExperimentType() {
		return null;
	}

	public String getExperimentDesign() {
		return null;
	}

	public String getAuthor() {
		return null;
	}

	public String getPublication() {
		return null;
	}

	public Array getArray(String name) {
		return null;
	}

	public Iterable<String> getArrayNames() {
		return null;
	}

}
