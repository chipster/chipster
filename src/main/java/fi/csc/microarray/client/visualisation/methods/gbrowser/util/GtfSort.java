package fi.csc.microarray.client.visualisation.methods.gbrowser.util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URISyntaxException;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.GtfLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.InMemoryIndex;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.LineDataSource;

/**
 * In memory sort, not for huge files.
 * 
 * @author klemela
 *
 */
public class GtfSort {
	
	public static void main(String[] args) throws FileNotFoundException, MalformedURLException, IOException, URISyntaxException, GBrowserException {
		
		//String fileString = System.getProperty("user.home") + "/chipster/Homo_sapiens.GRCh37.66.gtf";
		//String fileString = System.getProperty("user.home") + "/chipster/cufflinks-gtf/merged.gtf";
		String fileString = System.getProperty("user.home") + "/chipster/cufflinks-gtf/transcripts.gtf";
		
		File file = new File(fileString);				
		File outFile = new File(fileString.replace(".gtf", "-sort.gtf"));
		
		if (outFile.exists()) {
			System.err.println("Outpu file exists already!");
			System.exit(1);
		}
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
		
		InMemoryIndex index;
		
		DataUrl dataUrl = new DataUrl(file);

		index = new InMemoryIndex(new LineDataSource(dataUrl), new GtfLineParser());		
				
		for (String line : index.getFileLines()) {
			out.append(line);
			out.newLine();
			
			System.out.println(line.substring(0, 10));
		}
		
		out.flush();
		out.close();
	}
}