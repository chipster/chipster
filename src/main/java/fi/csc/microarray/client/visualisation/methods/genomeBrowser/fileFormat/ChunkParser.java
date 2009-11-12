package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher.ByteChunk;

public abstract class ChunkParser {

	protected ByteChunk chunk;
	protected FileDefinition fileDef;
	
	public ChunkParser(FileDefinition fileDef){
		
		this.fileDef = fileDef;
	}

	public void setChunk(ByteChunk chunk){
		this.chunk = chunk;
	}
	
	public abstract ChunkParser clone();

	public abstract long getReadCount();

	protected abstract long getPosition(long readIndex, Content field);

	protected abstract String getString(long pos);

	protected abstract long getLong(long pos);

	protected abstract float getFloat(long pos);
	
	public String getString(long readIndex, Content field){
		return getString(getPosition(readIndex, field));
	}

	public long getLong(long readIndex, Content field){
		return getLong(getPosition(readIndex, field));
	}

	public float getFloat(long readIndex, Content field){
		return getFloat(getPosition(readIndex, field));
	}
	
	public Object get(long pos, Type type){
		if( type == Type.FLOAT){
			return getFloat(pos);
		} else if( type == Type.LONG){
			return getLong(pos);
		} else if( type == Type.STRING){
			return getString(pos);
		} else {
			return null;
		}
	}
	
	public List<Object> getRead( long readIndex){
		
		List<Object> fields = new LinkedList<Object>();
		
		for(DataFieldDef field : fileDef){
			fields.add(this.get(getPosition(readIndex, field.content), 
					fileDef.getFieldDef(field.content).type));		
		}
		
		return fields;
	}

	public Map<Content, Object> getValues(long l,
			Collection<Content> requestedContents) {
		
		Map<Content, Object> values = new HashMap<Content, Object>();
		
		for(Content requestedContent : requestedContents){
			if(fileDef.getFieldDef(requestedContent) != null){
			values.put(requestedContent, 
					this.get(getPosition(l, requestedContent), 
							fileDef.getFieldDef(requestedContent).type));
			}
		}
		
		return values;
	}
}

