package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.Chromosome;

public abstract class ConstantRowLengthParser extends FileParser{

	private FileDefinition fileDef;
	private Integer rowLength;
	
	public ConstantRowLengthParser(FileDefinition fileDef){
		this.fileDef = fileDef;
	}
	
	public int getRowByteLength() {

		if(rowLength == null){
			rowLength = 0;
			for(ColumnDefinition dataFieldDef: fileDef){
				rowLength += dataFieldDef.length;
			}
		}
		return rowLength;
	}
	
	public Map<ColumnType, Object> getValues(long l,
			Collection<ColumnType> requestedContents) {

		Map<ColumnType, Object> values = new HashMap<ColumnType, Object>();

		for(ColumnType requestedContent : requestedContents){
			if(fileDef.getFieldDef(requestedContent) != null){
				values.put(requestedContent, 
						this.get(l, requestedContent));
			}
		}

		return values;
	}
	
	@Override
	public Object get(long rowIndex, ColumnType col) {
		
		ColumnDefinition fieldDef = fileDef.getFieldDef(col);

		byte[] byteValue = new byte[fieldDef.length];

		System.arraycopy(chunk.byteContent, 
				(int)((rowIndex - chunk.rowIndex) * getRowByteLength() + fieldDef.offset), 
				byteValue, 0, byteValue.length);
		
		String string = new String(byteValue).trim();

		if(col == ColumnType.STRAND) {
			return string.equalsIgnoreCase("r") || string.equals("-") ? 
					Strand.REVERSED : Strand.FORWARD;
			
		} else if(col == ColumnType.CHROMOSOME) {
			return new Chromosome(string.replace("chr", ""));

		} else if(fieldDef.type == Type.STRING) {
			return string;

		} else if(fieldDef.type == Type.FLOAT){
			return new Float(string);

		} else if(fieldDef.type == Type.LONG) {
			return new Long(string);
		}
		return null;
	}

	@Override
	public long getChunkRowCount() {
		return chunk.byteLength / getRowByteLength();
	}
}
