package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Generic parser for tab-separated value files.
 *
 */
public abstract class TsvParser extends FileParser {

		private FileDefinition fileDef;

		public TsvParser(FileDefinition fileDef) {
			this.fileDef = fileDef;
		}
		
		public String[] getFirstRow(Chunk chunk) {
			return chunk.getContent().substring(0, chunk.getContent().indexOf("\n")).split("\t");
		}
		
		public String[] getLastRow(Chunk chunk) {
			
			//minus two to convert from length to index and skip the last line change
			int lineStartIndex = chunk.getContent().lastIndexOf("\n", chunk.getContent().length() - 2);
			
			if (lineStartIndex < 0) {
				lineStartIndex = 0;
			}
			
			return chunk.getContent().substring(lineStartIndex+1, chunk.getContent().length()).split("\t");
		}
		
		@Override
		public BpCoordRegion getBpRegion(Chunk chunk) {
		    String[] firstRow = getFirstRow(chunk);
		    String[] lastRow = getLastRow(chunk);
			Long start = (Long)get(firstRow, ColumnType.BP_START);
			Long end;
            try {
                // Check if end is given
                end = (Long)get(lastRow, ColumnType.BP_END);
            } catch (RuntimeException e) {
                // Calculate from sequence
                end = (Long)get(lastRow, ColumnType.BP_START) +
                      (long)((String)get(lastRow, ColumnType.SEQUENCE)).length();
            }
			Chromosome startChr = (Chromosome)get(firstRow, ColumnType.CHROMOSOME);
			Chromosome endChr = (Chromosome)get(lastRow, ColumnType.CHROMOSOME);
			
			return new BpCoordRegion(start, startChr, end, endChr);
		}
		
		/**
		 * Fetch values of given column type from given array of columns.
		 * 
		 * @param cols - an array of columns in a row (as Strings).
		 * @param col - column type.
		 * @return
		 */
		public Object get(String[] cols, ColumnType col) {

			try {

				if (cols.length <= 1) {
					return null;
				}

				String string = cols[getFileDefinition().indexOf(col)].trim();

				ColumnDefinition fieldDef = getFileDefinition().getFieldDef(col);

				if (col == ColumnType.STRAND) {
					return string.equals("2") || string.equalsIgnoreCase("r") 
					|| string.equals("-") ? Strand.REVERSED	: Strand.FORWARD;

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
				
			} catch (Exception e) {
				throw new RuntimeException("error parsing columns: " + Arrays.toString(cols) + " (looking for: " + col + ")", e);
			}
		}
		
	    /**
         * Fetch all columns of given column types that are in a given chunk.
         * 
         * @param chunk - string chunk from a file.
         * @param requestedContents - a collection of column types to be fetched.
         * @return
         */
		@Override
		public List<RegionContent> getAll(Chunk chunk, Collection<ColumnType> requestedContents) {

			List<RegionContent> rows = new LinkedList<RegionContent>();
			

			for (String row : chunk.getContent().split("\n")) {
				
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
		
		public void setFileDefinition(FileDefinition fileDef) {
			this.fileDef = fileDef;
		}

		public FileDefinition getFileDefinition() {
			return fileDef;
		}
		
		@Override
		public long getDefaulChunkLength() {
			return 2*1024;
		}

		public long getHeaderLength(File file) {
			return 0;
		}
}
