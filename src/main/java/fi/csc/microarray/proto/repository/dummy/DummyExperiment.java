package fi.csc.microarray.proto.repository.dummy;

import java.util.LinkedList;

import fi.csc.microarray.proto.repository.Array;
import fi.csc.microarray.proto.repository.ExperimentBase;

public class DummyExperiment extends ExperimentBase {

	public String getName() {
		return "Dummy experiment";
	}
	
	public String getDescription() {
		return "Dummy description of the experiment";
	}

	public Array getArray(String name) {
		DummyArray arr = new DummyArray("cdc10");
		return arr;
	}

	public Iterable<String> getArrayNames() {
		LinkedList<String> arrs = new LinkedList<String>();
		arrs.add("cdc10");
		return arrs;
	}

	public String getUniqueIdentifier() {
		return "dummy1";
	}
}
