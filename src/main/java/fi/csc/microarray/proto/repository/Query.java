package fi.csc.microarray.proto.repository;

import java.util.LinkedList;

public class Query implements Comparable<Query> {
	
	private static final String SUBQUERY_SEPARATOR = " AND ";
	private static final String STEP_SEPARATOR = "/";
	protected String query;
	protected Repository repository;
	
	public Query(String query, Repository repository) {
		this.query = query;
		this.repository = repository;
	}

	public boolean isValid() {
		try {
			if (hasSubQueries()) {
				for (Query subQuery : subQueries()) {
					if (!subQuery.isValid()) {
						return false;
					}
				}
			} else {
				for (String step : steps()) {
					if (step.length() < 1) {
						return false;
					}
				}
			}
		} catch (Exception e) {
			return false;
		}
		return true;
	}
	
	@Override
	public String toString() {
		return query;
	}
	
	public Iterable<Query> subQueries() {
		LinkedList<Query> subQueries = new LinkedList<Query>();
		for (String subQuery : query.split(SUBQUERY_SEPARATOR)) {
			subQueries.add(new Query(subQuery.trim(), this.repository));
		}
		return subQueries;
	}
	
	public boolean hasSubQueries() {
		return query.contains(SUBQUERY_SEPARATOR);
	}
	
	public String getLastStep() {
		return query.substring(query.lastIndexOf(STEP_SEPARATOR)+1).trim();
	}
	
	public String getSecondLastStep() {
		String[] steps = query.split(STEP_SEPARATOR);
		return steps[steps.length-2];
	}
	
	public Repository getRepository() {
		return repository;
	}
	
	public Iterable<String> steps() {
		assert(!hasSubQueries());
		LinkedList<String> steps = new LinkedList<String>();
		for (String subQuery : query.split(STEP_SEPARATOR)) {
			steps.add(subQuery.trim());
		}
		return steps;
	}
	
	@Override
	public boolean equals(Object o) {
		if (o instanceof Query) {
			if (this.query.equals(((Query) o).toString())) {
				return true;
			}
		}
		return false;
	}
	
	public int compareTo(Query q) {
		return this.query.compareTo(q.toString());
	}
}