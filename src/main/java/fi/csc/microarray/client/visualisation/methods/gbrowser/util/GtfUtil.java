package fi.csc.microarray.client.visualisation.methods.gbrowser.util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.GtfLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.InMemoryIndex;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.LineDataSource;

/**
 * FIXME Find out why 400 MB file is only 150 MB after sorting. 
 * 
 * @author klemela
 */
public class GtfUtil {

	public static List<Feature> loadFile(File file) {
		InMemoryIndex index;
		List<Feature> rows = new LinkedList<Feature>();
		
		try {
			DataUrl dataUrl = new DataUrl(file);
			
			index = new InMemoryIndex(new LineDataSource(dataUrl), new GtfLineParser());
			
			GtfLineParser parser = new GtfLineParser();
			
			for (String line : index.getFileLines()) {
				parser.setLine(line);
				rows.add(new Feature(parser.getRegion()));
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
