package fi.csc.microarray.proto.repository;

import java.util.HashMap;

import fi.csc.microarray.exception.MicroarrayException;


public abstract class RepositoryBase extends Repository {

	@Override
	public boolean requiresAuthentication() {
		return false;
	}


	@Override
	public boolean authenticate(String username, String password) {
		throw new UnsupportedOperationException("repository does not require authentication");
	}
	
	@Override
	public Iterable<Experiment> executeQuery(Query query) throws MicroarrayException {
		
		HashMap<String, Experiment> intersection = new HashMap<String, Experiment>();
		
		for (Query subQuery : query.subQueries()) {
			for (Experiment experiment : executeSubQuery(subQuery)) {
				intersection.put(experiment.getUniqueIdentifier(), experiment);
			}
		}
		
		return intersection.values();
	}
	
	public abstract Iterable<Experiment> executeSubQuery(Query query) throws MicroarrayException;

}
