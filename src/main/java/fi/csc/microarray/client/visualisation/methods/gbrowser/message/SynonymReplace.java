package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class SynonymReplace {

	public static class Synonym {
		private String searchFor;
		private String replaceWith;

		public Synonym(String searchFor, String replaceWith) {
			this.searchFor = searchFor;
			this.replaceWith = replaceWith;
		}

		public String apply(String name) {
			if (searchFor.equals(name)) {
				return replaceWith;
			}
			return name;
		}

		public String revoke(String name) {
			if (replaceWith.equals(name)) {
				return searchFor;
			}
			return name;
		}

		public String getReplaceWith() {
			return replaceWith;
		}

		public String getSearchFor() {
			return searchFor;
		}
	}

	public List<Synonym> list;

	public SynonymReplace(Synonym[] array) {
		this.list = Arrays.asList(array);
	}

	public SynonymReplace() {
		list = new LinkedList<>();
	}

	public String apply(String name) {
		for (Synonym synonym : list) {
			name = synonym.apply(name);
		}
		return name;
	}

	public String revoke(String name) {
		for (Synonym synonym : list) {
			name = synonym.revoke(name);
		}
		return name;
	}

	public List<Synonym> getSynonyms() {
		return list;
	}

	public void add(Synonym synonym) {
		list.add(synonym);		
	}
}