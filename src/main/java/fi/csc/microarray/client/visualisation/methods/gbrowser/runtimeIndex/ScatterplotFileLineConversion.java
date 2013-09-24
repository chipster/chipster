package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.Iterator;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;

/**
 * This class converts bed files to BedLine objects.
 *  
 * Results are in RegionContent objects. These objects contain:
 * <ul>
 * <li>the region
 * <li>an unique line identifier stored with key DataType.ID
 * <li>the BedLine object stored with key DataType.VALUE
 * </ul> 
 * 
 * @author klemela
 *
 */
public class ScatterplotFileLineConversion extends FileLineConversion {

	private Float minScatterplotValue;
	private Float maxScatterplotValue;

	public ScatterplotFileLineConversion(DataUrl data, AbstractTsvLineParser parser, GBrowser browser) throws URISyntaxException, IOException, GBrowserException {
		super(data, parser, browser);
		
		udpatesScatterplotValues();
	}

	public Float getMinScatterplotValue() {
		return minScatterplotValue;
	}

	public Float getMaxScatterplotValue() {
		return maxScatterplotValue;
	}
	
	private void udpatesScatterplotValues() throws IOException, GBrowserException {
		
		boolean notSet = true;
		minScatterplotValue = Float.MAX_VALUE;
		maxScatterplotValue = Float.NEGATIVE_INFINITY;
		
		Iterator<String> iter = super.getIndex().getFileLineIterator();
		
		while (iter.hasNext()) {
			LineParser parser = getParser();
			parser.setLine(iter.next());
			if (parser.isContentLine()) {
				ScatterplotValue line = (ScatterplotValue)parser.getFileLine();

				float value = line.getScatterplotValue();

				//Cufflinks produces scores that float can not present (e.g. -1.79769e+308). Do not let those to trample min/max values.  
				if (!Float.isNaN(value) && !Float.isInfinite(value)) {
					if (value > maxScatterplotValue && value < Float.MAX_VALUE) {
						maxScatterplotValue = value;
					}

					if (value < minScatterplotValue && value > -Float.MAX_VALUE) {
						minScatterplotValue = value;
					}
					notSet = false;
				}

			}
		}
		
		if (notSet) {
			minScatterplotValue = null;
			maxScatterplotValue = null;
		}		
	}
}
