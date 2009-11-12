package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import java.util.Collection;
import java.util.LinkedList;


public class FileDefinition extends LinkedList<DataFieldDef> {
	
	public FileDefinition(Collection<DataFieldDef> content){
		this.addAll(content);
	}
	
	public DataFieldDef getFieldDef(Content type){
		for(DataFieldDef field: this){
			if ( field.content == type){
				return field;
			}
		}
		
		return null;
	}
	
	public int indexOf(Content content){
		int i = 0;
		for(DataFieldDef field : this){
			if(field.content == content){
				return i;
			}
			i++;
		}
		return -1;
	}	
}
