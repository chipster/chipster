package fi.csc.microarray.gbrowser.index;


import java.io.File;
import java.io.IOException;

import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.GeneIndexParser;

public class GetData {
	
	private static final File URL_ROOT;

	static {
			URL_ROOT = new File("/home/zukauska/chipster-share/ngs/annotations");
	}

	public GetData (){
		
	}
	
	public static void main(String args[]) throws IOException{
		ChunkDataSource cds = new ChunkDataSource(URL_ROOT, "Homo_sapiens.NCBI36.54_genes.tsv", new GeneIndexParser());
		byte[] a; 
		int asdf = 1234;
		a = new byte[asdf];
		cds.read(0, a);
		
		System.out.println(cds));
		//System.out.println("asdf"+cds.read(0, ));
	}
}
