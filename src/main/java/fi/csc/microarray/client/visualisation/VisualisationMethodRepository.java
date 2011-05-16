package fi.csc.microarray.client.visualisation;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;

public class VisualisationMethodRepository {
	
	private LinkedList<VisualisationMethod> visualisationMethods = new LinkedList<VisualisationMethod>(); 

	
	public void addVisualisationMethods(VisualisationMethod[] newMethods) {
		this.visualisationMethods.addAll(Arrays.asList(newMethods));
	}
	
	public List<VisualisationMethod> getVisualisationMethods() {
		return visualisationMethods;
	}
	
	public VisualisationMethod getDefaultVisualisationFor(DataBean dataBean) throws IOException, MicroarrayException {
		for (VisualisationMethod method : getOrderedDefaultCandidates()) {
			if (method != VisualisationMethod.NONE && method.isApplicableTo(dataBean)) {
				return method;
			}
		}
		return null;
	}
	
	public Iterable<VisualisationMethod> getOrderedDefaultCandidates() {
		
		LinkedList<VisualisationMethod> orderedDefaultCandidates = new LinkedList<VisualisationMethod>();
		for (VisualisationMethod method : visualisationMethods) {
			if (method.getOrderNumber() >= 0) {
				orderedDefaultCandidates.add(method);
			}
		}
		Collections.sort(orderedDefaultCandidates, new Comparator<VisualisationMethod>() {
			public int compare(VisualisationMethod method1, VisualisationMethod method2) {
				return new Integer(method2.getOrderNumber()).compareTo(new Integer(method1.getOrderNumber()));
			}
		});

		return orderedDefaultCandidates;
	}	

}
