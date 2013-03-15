package fi.csc.microarray.client.visualisation.methods.gbrowser.util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.LineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.stack.GtfToFeatureConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.stack.InMemoryIndex;
import fi.csc.microarray.client.visualisation.methods.gbrowser.stack.GtfLineParser;

/**
 * FIXME Find out why 400 MB file is only 150 MB after sorting. 
 * 
 * 
 * @author klemela
 */
public class GtfUtil {

	public static List<RegionContent> loadFile(File file) {
		InMemoryIndex index;
		List<RegionContent> rows = new LinkedList<RegionContent>();
		
		try {
			index = new InMemoryIndex(new LineDataSource(file.toURI().toURL(), GtfToFeatureConversion.class), new GtfLineParser());
			
			GtfLineParser parser = new GtfLineParser();
			
			for (String line : index.getFileLines()) {
				parser.setLine(line);
				rows.add(new RegionContent(parser.getRegion()));
			}
			
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (MalformedURLException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (URISyntaxException e) {
			e.printStackTrace();
		}
		
		return rows;
	}
}
