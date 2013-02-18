package fi.csc.microarray.client.visualisation.methods.gbrowser.util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URISyntaxException;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.LineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.stack.InMemoryIndex;
import fi.csc.microarray.client.visualisation.methods.gbrowser.stack.StackGtfParser;

/**
 * In memory sort, not for huge files.
 * 
 * @author klemela
 *
 */
public class GtfSort {
	
	public static void main(String[] args) throws FileNotFoundException, MalformedURLException, IOException, URISyntaxException, GBrowserException {
		
		File file = new File(System.getProperty("user.home") + "/chipster/Homo_sapiens.GRCh37.66.gtf");
		
		File outFile = new File(System.getProperty("user.home") + "/chipster/Homo_sapiens.GRCh37.66-sort.gtf");
		BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
		
		InMemoryIndex index;

		index = new InMemoryIndex(new LineDataSource(file.toURI().toURL(), null), new StackGtfParser());		
				
		for (String line : index.getFileLines()) {
			out.append(line);
			out.newLine();
			
			System.out.println(line.substring(0, 10));
		}
		
		out.flush();
		out.close();
	}
}