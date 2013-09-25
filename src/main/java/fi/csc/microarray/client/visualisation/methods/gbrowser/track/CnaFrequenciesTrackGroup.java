package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.CnaConversion;

public class CnaFrequenciesTrackGroup extends CnaTrackGroup {

	public CnaFrequenciesTrackGroup(GBrowserView view,
			CnaConversion dataThread, LinkedList<String> sampleNames, String title) {
		super(view, dataThread, sampleNames, title + " frequencies", false);
	}
	
	@Override
	public void addTracks() {			
		super.addTracks();
		
		if (!isMinimized()) {
			TitleTrack title2 = new TitleTrack("loss frequency", Color.black, GBrowserConstants.SCATTERPLOT_TITLE_COLOR);
			title2.setView(view);
			addTrack(title2);

			ScatterplotTrack lossFreq = new ScatterplotTrack(GBrowserConstants.BED_COLOR, 100, 0f, 1.0f, DataType.LOSS);
			lossFreq.setView(view);
			lossFreq.addDataThread(dataThread);
			addTrack(lossFreq);

			TitleTrack title3 = new TitleTrack("gain frequency", Color.black, GBrowserConstants.SCATTERPLOT_TITLE_COLOR);
			title3.setView(view);
			addTrack(title3);

			ScatterplotTrack gainFreq = new ScatterplotTrack(GBrowserConstants.BED_COLOR, 100, 0f, 1.0f, DataType.GAIN);
			gainFreq.setView(view);
			gainFreq.addDataThread(dataThread);
			addTrack(gainFreq);
		}
	}
}
