package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.net.URL;

public class GenomeInfo {
	
	private String species;
	private String version;
	private URL ensemblBrowserUrl;
	private URL ucscBrowserUrl;
	private String sortId;
	
	public String getSpecies() {
		return species;
	}
	
	public void setSpecies(String species) {
		this.species = species;
	}
	
	public String getVersion() {
		return version;
	}
	
	public void setVersion(String version) {
		this.version = version;
	}	
	
	public URL getEnsembl() {
		return ensemblBrowserUrl;
	}
	
	public void setEnsemblBrowserUrl(URL ensemblBrowserUrl) {
		this.ensemblBrowserUrl = ensemblBrowserUrl;
	}
	
	public URL getBrowserUrl() {
		return ucscBrowserUrl;
	}
	
	public void setUcscBrowserUrl(URL ucscBrowserUrl) {
		this.ucscBrowserUrl = ucscBrowserUrl;
	}

	public String getSortId() {
		return sortId;
	}

	public void setSortId(String sortId) {
		this.sortId = sortId;
	}
}