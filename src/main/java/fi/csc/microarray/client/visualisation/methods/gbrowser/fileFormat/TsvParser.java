package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public abstract class TsvParser extends FileParser {

		private FileDefinition fileDef;

		public TsvParser(FileDefinition fileDef) {
			this.fileDef = fileDef;
		}
		
		public String[] getFirstRow() {
			return chunk.substring(0, chunk.indexOf("\n")).split("\t");
		}
		
		public String[] getLastRow() {
			
			//minus two to convert from length to index and skip the last line change
			return chunk.substring(chunk.lastIndexOf("\n", chunk.length() - 2), chunk.length()).split("\t");
		}
		
		@Override
		public BpCoordRegion getBpRegion() {
			Long start = (Long)get(getFirstRow(), ColumnType.BP_START);
			Long end = (Long)get(getLastRow(), ColumnType.BP_START);
			Chromosome chr = (Chromosome)get(getFirstRow(), ColumnType.CHROMOSOME);
			
			return new BpCoordRegion(start, end, chr);
		}
		
		
		public Object get(String[] cols, ColumnType col) {
			
			
			if(cols.length <= 1) {
				return null;
			}
			
			String string = cols[fileDef.indexOf(col)].trim();

			ColumnDefinition fieldDef = fileDef.getFieldDef(col);
			
			if (col == ColumnType.STRAND) {
				return string.equalsIgnoreCase("r") || string.equals("-") ? Strand.REVERSED
						: Strand.FORWARD;

			} else if (col == ColumnType.CHROMOSOME) {
				return new Chromosome(string.replace("chr", ""));

			} else if (fieldDef.type == Type.STRING) {
				return string;

			} else if (fieldDef.type == Type.FLOAT) {
				return new Float(string);

			} else if (fieldDef.type == Type.LONG) {

				if (string.length() > 0) {
					return new Long(string);
				} else {
					return Long.MIN_VALUE;
				}
			}
			return null;
		}

		@Override
		public List<RegionContent> getAll(Collection<ColumnType> requestedContents) {

			List<RegionContent> rows = new LinkedList<RegionContent>();
			

			for (String row : chunk.split("\n")) {
				
				Map<ColumnType, Object> values = new HashMap<ColumnType, Object>();
				
				String[] cols = row.split("\t");
				
				for (ColumnType requestedContent : requestedContents) {
							
					values.put(requestedContent, this.get(cols, requestedContent));					
				}
				
				Long start = (Long)get(cols, ColumnType.BP_START);
				Long end = (Long)get(cols, ColumnType.BP_END);
				Chromosome chr = (Chromosome)get(cols, ColumnType.CHROMOSOME);
		
				rows.add(new RegionContent(new BpCoordRegion(start, end, chr), values));

			}
			
			return rows;
		}

		public FileDefinition getFileDefinition() {
			return fileDef;
		}
		
		@Override
		public long getDefaulChunkLength() {

			return 8*1024;
		}
}
