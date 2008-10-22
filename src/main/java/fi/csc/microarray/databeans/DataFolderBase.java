package fi.csc.microarray.databeans;

import fi.csc.microarray.util.Strings;


public abstract class DataFolderBase implements DataFolder {

	public DataFolder getChildFolder(String name) {
		for (DataItem child : getChildren()) {
			if (child instanceof DataFolder && child.getName().equals(name)) {
				return (DataFolder)child;
			}
		}
		return null;
	}
	
	public String toString() {
		return getName();
	}

	public String toStringRecursively(int i) {

		String s = (Strings.repeat("  ", i) + "[" + getName() + "]");
		
		for (DataItem child : getChildren()) {
			s += ("\n" + child.toStringRecursively(i+1));
		}
		
		return s;
	}
}
