package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;

import javax.swing.SwingUtilities;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.CnaRow;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.IndexKey;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.CnaRow.Sample;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.UnsortedDataException;

/**
 * This conversion class Parses tsv files with CnaLineParser, creates CnaRow objects and packages those as
 * RegionContent objects.  
 * 
 * @author klemela
 *
 */
public class CnaConversion extends SingleThreadAreaRequestHandler {

	private Index index;

	private CnaLineParser parser;

	public CnaConversion(DataSource file, final GBrowser browser) {
	    
		super(null, null);

		this.parser = new CnaLineParser();
		
		try {
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
	protected void processAreaRequest(AreaRequest request) {
						
		super.processAreaRequest(request);	
		
		if (request.getStatus().poison) {
			return;
		}
		
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
						
			lines = index.getFileLines(new AreaRequest(requestRegion, request.getRequestedContents(), request.getStatus()));
			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (GBrowserException e) {
			e.printStackTrace();
		}
		
		List<RegionContent> list = new LinkedList<RegionContent>();
		
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
			
			LinkedHashMap<ColumnType, Object> valueMap = new LinkedHashMap<ColumnType, Object>();			
			
			valueMap.put(ColumnType.ID, id);
			valueMap.put(ColumnType.VALUE, row);
			valueMap.put(ColumnType.LOSS, row.getLossFreg());
			valueMap.put(ColumnType.GAIN, row.getGainFreg());
			
			//Add logRatios in general format to make it possible to view them with ScatterploTrack			
			valueMap.put(ColumnType.FLOAT_LIST, logRatioValues);
			
			RegionContent regionContent = new RegionContent(region, valueMap);
			
			list.add(regionContent);
		}	
		
		super.createAreaResult(new AreaResult(request.getStatus(), list));
	}
}
