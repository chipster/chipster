package fi.csc.microarray.proto.repository.schema;

import java.util.ArrayList;
import java.util.List;

public class ParameterInstance extends ParameterItem {
	
	private List<String> synonyms;
	
	public ParameterInstance(String name) {
		super(name);
		this.synonyms = null;
	}
	
	public void addSynonym(String synonym) {
		if (synonym != null) {
			if (synonyms == null) {
				synonyms = new ArrayList<String>();
			}
			synonyms.add(synonym);
		}
	}
	
	public String[] getSynonyms() {
		if (synonyms != null) {
			return (String[]) synonyms.toArray(new String[] { });
		} else {
			return null;
		}
	}
	
	public int getSynonymCount() {
		if (synonyms != null) {
			return synonyms.size();
		} else {
			return 0;
		}
	}
	
	public boolean hasSynonyms() {
		return (synonyms != null && synonyms.size() > 0);
	}
}