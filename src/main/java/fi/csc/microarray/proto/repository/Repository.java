package fi.csc.microarray.proto.repository;

import java.util.LinkedList;
import java.util.Properties;

import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.proto.repository.schema.ParameterClass;

public abstract class Repository {
	
	
	//
	// static implementation part
	//
	
	private static final String REPOSITORY_PROPERTY_PREFIX = "nami.repository.available.";

	public static Iterable<Repository> getAvailableRepositories() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		LinkedList<Repository> repositories = new LinkedList<Repository>();
		Properties properties = System.getProperties();
		for (Object name : properties.keySet()) {
			if (name instanceof String && ((String)name).startsWith(REPOSITORY_PROPERTY_PREFIX)) {
				String className = (String)properties.get(name);
				Repository repository = (Repository)Class.forName(className).newInstance();
				repositories.add(repository);
			}
		}
		return repositories;
	}
	
	
	//
	// non-static interface part
	//

	public abstract boolean requiresAuthentication();
	
	public abstract boolean authenticate(String username, String password);

	/**
	 * Identifier for the actual implementation behind the interface.
	 */
	public abstract String getType();

	/**
	 * Identifier (e.g. URL) for the physical data source.
	 */
	public abstract String getIdentifier();
	
	public abstract ParameterClass getRootClass() throws MicroarrayException;
	
	public abstract int getClassHierarchyDepth();
	
	public abstract Iterable<Experiment> executeQuery(Query query) throws MicroarrayException;

	public abstract String toString();
}
