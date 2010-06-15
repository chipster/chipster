package fi.csc.microarray.databeans.features;

import java.util.LinkedList;

import fi.csc.microarray.databeans.Dataset;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.util.LookaheadStringReader;

public class RequestExecuter {
	
	private DataManager manager;

	public RequestExecuter(DataManager manager) {
		this.manager = manager;
	}
	
	public Feature execute(String request, Dataset data) {
		
		if (request.startsWith("/")) {
			return executeFeature(request, data);
		} else {
			return executeModifier(request, data);
		}
	}

	private Feature executeModifier(String request, Dataset data) {
		LookaheadStringReader requestReader = new LookaheadStringReader(request);
		
		String modifierName = requestReader.readTo("(");
		Modifier modifier = manager.fetchModifier(modifierName);
		
		requestReader.read(); // read '('
		
		String modifiedExpression = requestReader.readToLast(")");
		
		Feature feature = execute(modifiedExpression, data);
		LinkedList<Feature> fl = new LinkedList<Feature>();
		fl.add(feature);
		modifier.setInputs(fl);
		
		requestReader.read(); // read ')'
		
		if (!requestReader.isAtEnd()) {
			throw new IllegalArgumentException("malformed request: " + request);
		}
		
		return modifier.getOutput();
	}

	private Feature executeFeature(String request, Dataset data) {
		return manager.fetchFeature(request, data);
	}
}
