package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserPlot;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.stack.CnaConversion;

/**
 * Utility class for creating predefined {@link TrackGroup} objects.  
 *
 *@author Petri Klemelä, Aleksi Kallio
 *
 */
public class TrackFactory {
	
    public static TrackGroup getGeneTrackGroup(GBrowserPlot genomePlot, AreaRequestHandler annotationDataSource, AreaRequestHandler repeatDataSource, boolean isUserData) throws FileNotFoundException {
        
		GBrowserView dataView = genomePlot.getDataView();
		
		TrackGroup geneGroup = new GeneTrackGroup(dataView, annotationDataSource, repeatDataSource, isUserData);
			    
	    return geneGroup;
	}
    
	public static TrackGroup getThinSeparatorTrackGroup(GBrowserPlot genomePlot) {
		GBrowserView view = genomePlot.getDataView();
		TrackGroup group = new TrackGroup(view);
		SeparatorTrack separator = new SeparatorTrack(Color.LIGHT_GRAY, 3, 0, Long.MAX_VALUE); 
		separator.setView(view);
		group.addTrack(separator);
		return group;
	}

	public static TrackGroup getThickSeparatorTrackGroup(GBrowserPlot genomePlot) {
		GBrowserView view = genomePlot.getDataView();
		TrackGroup group = new TrackGroup(view);
		
		SeparatorTrack separator1 = new SeparatorTrack3D(0, Long.MAX_VALUE, false);		
		separator1.setView(view);	
		group.addTrack(separator1);
		
		EmptyTrack empty = new EmptyTrack(2);
		empty.setView(view);		
		group.addTrack(empty);
		
		SeparatorTrack separator2 = new SeparatorTrack3D(0, Long.MAX_VALUE, true);		
		separator2.setView(view);	
		group.addTrack(separator2);
		
		return group;
	}
	
	public static TrackGroup getReadTrackGroup(GBrowserPlot genomePlot, AreaRequestHandler details, AreaRequestHandler estimate, 
			AreaRequestHandler seqFile, String title)	
	        throws FileNotFoundException, MalformedURLException {
	
		GBrowserView dataView = genomePlot.getDataView();
		
		// Group containing tracks for this data source
		ReadTrackGroup readGroup = new ReadTrackGroup(dataView, details, estimate, seqFile, title);
		readGroup.initialise();
        
        return readGroup;
	}

	public static TrackGroup getReadSummaryTrackGroup(GBrowserPlot genomePlot,  AreaRequestHandler details, AreaRequestHandler estimate, 
			AreaRequestHandler seqFile, String title, AreaRequestHandler summaryDataSource)
	        throws FileNotFoundException, MalformedURLException {
	
		GBrowserView dataView = genomePlot.getDataView();
		
		// Group containing tracks for this data source
		ReadSummaryTrackGroup readGroup = new ReadSummaryTrackGroup(dataView, details, estimate, seqFile, title, summaryDataSource);
		readGroup.initialise();
        
        return readGroup;
	}
	
	public static TrackGroup getPeakTrackGroup(GBrowserPlot plot, AreaRequestHandler areaRequestHandler) {
		GBrowserView dataView = plot.getDataView();

		
		PeakTrack annotation = new PeakTrack(GBrowserConstants.BED_COLOR, 0, Long.MAX_VALUE);
		annotation.setView(dataView);
		annotation.addAreaRequestHandler(areaRequestHandler);
		
		return new TrackGroup(annotation);
	}

	public static TrackGroup getCytobandTrackGroup(GBrowserPlot plot, AreaRequestHandler cytobandData) {
		
		CytobandTrack overviewCytobands = new CytobandTrack(false);
		overviewCytobands.setView(plot.getOverviewView());
		overviewCytobands.addAreaRequestHandler(cytobandData);
		
		return new TrackGroup(overviewCytobands);
	}

	public static TitleTrack getTitleTrack(GBrowserPlot genomePlot, String title) {
		GBrowserView dataView = genomePlot.getDataView();
		
		TitleTrack titleTrack = new TitleTrack(title, Color.black);
		titleTrack.setView(dataView);
		return titleTrack;
	}

	public static TrackGroup getCnaTrackGroup(GBrowserPlot plot,
			CnaConversion conversion, LinkedList<String> sampleNames, boolean showFrequencies, boolean showCalls, boolean showLogratios) {
		
		GBrowserView view = plot.getDataView();
		TrackGroup group = new TrackGroup(view);
		
		if (showFrequencies) {
			TitleTrack title2 = new TitleTrack("loss frequency", Color.black, GBrowserConstants.SCATTERPLOT_TITLE_COLOR);
			title2.setView(view);
			group.addTrack(title2);

			ScatterplotTrack lossFreq = new ScatterplotTrack(GBrowserConstants.BED_COLOR, 100, 0f, 1.0f, ColumnType.LOSS, 0, Long.MAX_VALUE);
			lossFreq.setView(view);
			lossFreq.addAreaRequestHandler(conversion);
			group.addTrack(lossFreq);

			TitleTrack title3 = new TitleTrack("gain frequency", Color.black, GBrowserConstants.SCATTERPLOT_TITLE_COLOR);
			title3.setView(view);
			group.addTrack(title3);
			
			ScatterplotTrack gainFreq = new ScatterplotTrack(GBrowserConstants.BED_COLOR, 100, 0f, 1.0f, ColumnType.GAIN, 0, Long.MAX_VALUE);
			gainFreq.setView(view);
			gainFreq.addAreaRequestHandler(conversion);
			group.addTrack(gainFreq);
		}
				
		for (int i = 0; i < sampleNames.size(); i++) {
			
			String name = sampleNames.get(i);

			if (showCalls) {
				
				TitleTrack title = new TitleTrack(name, Color.black);
				title.setView(view);
				group.addTrack(title);
				
				CnaFlagTrack flag = new CnaFlagTrack(GBrowserConstants.BED_COLOR, i, Color.RED, 0, Long.MAX_VALUE);
				flag.setView(view);
				flag.addAreaRequestHandler(conversion);
				group.addTrack(flag);

				SeparatorTrack separator1 = new SeparatorTrack(Color.gray, 1, 0, Long.MAX_VALUE);
				separator1.setView(view);
				group.addTrack(separator1);
			}
			
			if (showLogratios) {
				
				TitleTrack title = new TitleTrack(name, Color.black, GBrowserConstants.SCATTERPLOT_TITLE_COLOR);
				title.setView(view);
				group.addTrack(title);			
				
				ScatterplotTrack logRatio = new ScatterplotTrack(GBrowserConstants.BED_COLOR, 24, -1.0f, 1.0f, i, 0, Long.MAX_VALUE);				
				logRatio.setView(view);
				logRatio.addAreaRequestHandler(conversion);		
				group.addTrack(logRatio);
			}
		}		

		return group;
	}
}
