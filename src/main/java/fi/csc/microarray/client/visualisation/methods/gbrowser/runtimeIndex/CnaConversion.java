package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;

import javax.swing.SwingUtilities;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.CnaRow;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.IndexKey;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.CnaRow.Sample;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.UnsortedDataException;

/**
 * This conversion class Parses tsv files with CnaLineParser, creates CnaRow objects and packages those as
 * Feature objects.  
 * 
 * @author klemela
 *
 */
public class CnaConversion extends DataThread {

	private Index index;

	private CnaLineParser parser;

	private LinkedList<String> sampleNames;

	public CnaConversion(RandomAccessLineDataSource file, final GBrowser browser) {
			    
		super(browser, file);
		
		try {
			//read sample names
			file.setLineReaderPosition(0);
			String header = file.getNextLine();
			parser = new CnaLineParser();
			parser.setLine(header);
			
			this.sampleNames = parser.getSampleNames();
			
			//create index
			this.index = new BinarySearchIndex(file, parser);
		
		} catch (final UnsortedDataException e) {
			SwingUtilities.invokeLater(new Runnable() {
				
				@Override
				public void run() {				
					browser.showDialog("Unsorted data", e.getMessage(), null, true, false, true, true);
				}
			});
			index = null;
		} catch (IOException e) {
			e.printStackTrace();
		} catch (GBrowserException e) {
			e.printStackTrace();
		}
	}

	@Override
	protected void processDataRequest(DataRequest request) {						
		
		if (index == null) {
			return;
		}
		
		long start = request.start.bp;
		long end = request.end.bp;
		
		//TODO Be more clever to avoid getting so much useless data, get now the whole beginning of the chromosome, because cna regions can be tens of megabytes long
		start = 1;
		
		Region requestRegion = new Region(start, end, request.start.chr);		
		
		TreeMap<IndexKey, String> lines = null;
		try {		
						
			lines = index.getFileLines(new DataRequest(requestRegion, request.getRequestedContents(), request.getStatus()));
			
		} catch (IllegalArgumentException e) {
			e.printStackTrace();
			return;
		} catch (IOException e) {
			e.printStackTrace();
		} catch (GBrowserException e) {
			e.printStackTrace();
		}
		
		List<Feature> list = new LinkedList<Feature>();
		
		for (Entry<IndexKey, String> lineEntry : lines.entrySet()) {

			parser.setLine(lineEntry.getValue());					
			
			Region region = parser.getRegion();
			
			if (!region.intersects(request)) {
				continue;
			}
			
			CnaRow row = new CnaRow();
			row.setRegion(region);
			
			Float gainFreq = parser.getGainFreq();
			if (gainFreq != null) {
				row.setGainFreg(gainFreq);
			}
			
			Float lossFreq = parser.getLossFreq();
			if (lossFreq != null) {
				row.setLossFreg(parser.getLossFreq());
			}
			
			List<String> sampleNames = parser.getSampleNames();
			List<Float> flagValues = parser.getFlagValues();
			List<Float> logRatioValues = parser.getLogRatioValues();
			
			LinkedList<Sample> samples = new LinkedList<Sample>();
			
			for (int i = 0; i < sampleNames.size(); i++) {
				
				Sample sample = new Sample();
				
				sample.setName(sampleNames.get(i));
				
				if (flagValues.size() > 0) {
					sample.setFlag(flagValues.get(i));
				}
				
				sample.setLogRatio(logRatioValues.get(i));
				
				samples.add(sample);
			}
			
			row.setSamples(samples);
			
			IndexKey id = lineEntry.getKey();
			
			LinkedHashMap<DataType, Object> valueMap = new LinkedHashMap<DataType, Object>();			
			
			valueMap.put(DataType.ID, id);
			valueMap.put(DataType.VALUE, row);
			valueMap.put(DataType.LOSS, row.getLossFreg());
			valueMap.put(DataType.GAIN, row.getGainFreg());
			
			//Add logRatios in general format to make it possible to view them with ScatterploTrack			
			valueMap.put(DataType.FLOAT_LIST, logRatioValues);
			
			Feature regionContent = new Feature(region, valueMap);
			
			list.add(regionContent);
		}	
		
		super.createDataResult(new DataResult(request.getStatus(), list));
	}
	
	public LinkedList<String> getSampleNames() {
		return sampleNames;
	}
}
