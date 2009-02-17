package fi.csc.microarray.analyser.ws;

import org.w3c.dom.Document;

import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;
import fi.csc.microarray.analyser.ws.resultgraph.ResultGraph;

public abstract class SoapAnalysisJob extends JavaAnalysisJobBase {

	public abstract Document getEnvelope(Iterable<String> query);

	public abstract Document getResultTable(ResultGraph resultGraph);
}
