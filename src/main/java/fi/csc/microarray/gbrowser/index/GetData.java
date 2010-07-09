package fi.csc.microarray.gbrowser.index;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.GeneParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class GetData {
	
	private static final File URL_ROOT;

	static {
			URL_ROOT = new File("/home/zukauska/chipster-share/ngs/annotations");
	}
	

	public GetData (){
		
	}
	
	public List<RegionContent> read(){

        BufferedReader bufferedReader = null;
        ArrayList<GeneIndexDataType> indexList = null;
        
        try {
        	File geneFile = new File(URL_ROOT+"/Homo_sapiens.NCBI36.54_genes.tsv");
			ChunkDataSource data = new ChunkDataSource(
					geneFile,
					new GeneParser());
			byte[] fileChunk = new byte[(int)geneFile.length()];
			data.read(0, fileChunk);
			List<ColumnType> columns = Arrays.asList(new ColumnType[] {
					ColumnType.CHROMOSOME,
					ColumnType.BP_START,
					ColumnType.BP_END,
					ColumnType.DESCRIPTION});
			
			return  data.getFileParser().getAll(new String(fileChunk), columns);
			
//			for (RegionContent gene : genes) {
//				 System.out.println(gene.values.get(ColumnType.DESCRIPTION));
//				 System.out.println(gene.region.start);
//				 break;
//			}
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return new ArrayList();
        
//        try {
//            
//            bufferedReader = new BufferedReader(new FileReader(URL_ROOT+"/Homo_sapiens.NCBI36.54_genes.tsv"));
//            indexList = new ArrayList<GeneIndexDataType>();
//            
//            String line = null;
//            while (((line = bufferedReader.readLine()) != null)) {
//            	String[] split = line.split("	");
//            	indexList.add(new GeneIndexDataType(Long.parseLong(split[0]), Long.parseLong(split[1]), Long.parseLong(split[2]), split[4]));
//            }
//            
//        } catch (FileNotFoundException ex) {
//            ex.printStackTrace();
//        } catch (IOException ex) {
//            ex.printStackTrace();
//        } finally {
//            try {
//                if (bufferedReader != null)
//                    bufferedReader.close();
//            } catch (IOException ex) {
//                ex.printStackTrace();
//            }
//        }
//		return indexList;
	}
}
