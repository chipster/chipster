package fi.csc.chipster.web.list;

import java.util.HashMap;
// import com.google.gwt.dev.util.collect.HashMap;

public enum Type {

	GENERIC (0, "GENERIC"),
	GENE_EXPRS (1, "GENE_EXPRS");

	static HashMap<String, Integer> map = new HashMap<String, Integer>();
	static {
		for(int i = 0; i < Type.values().length; i++) {
			Type type = Type.values()[i];
			map.put(type.name, type.id);
		}
	}
	
	private int id;
	private String name;
	
	private Type(int id, String name) {
		this.id = id;
		this.name = name;
	}
	
	public String getName() {
		return name;
	}
	
	public int getId() {
		return id;
	}
	
	public int getIdByName(String name) {
		return map.get(name.toUpperCase());
	}
}
