package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
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
	
	public static TrackGroup getReadTrackGroup(GBrowserPlot genomePlot, AreaRequestHandler userData, AreaRequestHandler seqFile, String title)
	        throws FileNotFoundException, MalformedURLException {
	
		GBrowserView dataView = genomePlot.getDataView();
		
		// Group containing tracks for this data source
		ReadTrackGroup readGroup = new ReadTrackGroup(dataView, userData, seqFile, title);
		readGroup.initialise();
        
        return readGroup;
	}

	public static TrackGroup getReadSummaryTrackGroup(GBrowserPlot genomePlot, AreaRequestHandler userData,
			AreaRequestHandler seqFile, String title, AreaRequestHandler summaryDataSource)
	        throws FileNotFoundException, MalformedURLException {
	
		GBrowserView dataView = genomePlot.getDataView();
		
		// Group containing tracks for this data source
		ReadSummaryTrackGroup readGroup = new ReadSummaryTrackGroup(dataView, userData, seqFile, title, summaryDataSource);
		readGroup.initialise();
        
        return readGroup;
	}
	
	public static TrackGroup getPeakTrackGroup(GBrowserPlot plot, AreaRequestHandler areaRequestHandler) {
		GBrowserView dataView = plot.getDataView();

		
		PeakTrack annotation = new PeakTrack(GBrowserConstants.BED_COLOR, 0, Long.MAX_VALUE);
		annotation.setView(dataView);
		annotation.setAreaRequestHandler(areaRequestHandler);
		
		return new TrackGroup(annotation);
	}

	public static TrackGroup getCytobandTrackGroup(GBrowserPlot plot, AreaRequestHandler cytobandData) {
		
		CytobandTrack overviewCytobands = new CytobandTrack(false);
		overviewCytobands.setView(plot.getOverviewView());
		overviewCytobands.setAreaRequestHandler(cytobandData);
		
		return new TrackGroup(overviewCytobands);
	}

	public static TitleTrack getTitleTrack(GBrowserPlot genomePlot, String title) {
		GBrowserView dataView = genomePlot.getDataView();
		
		TitleTrack titleTrack = new TitleTrack(title, Color.black);
		titleTrack.setView(dataView);
		return titleTrack;
	}

	public static TrackGroup getCnaTrackGroup(GBrowserPlot plot,
			CnaConversion conversion) {
		
		SeparatorTrack separator1 = new SeparatorTrack(Color.gray, 1, 0, Long.MAX_VALUE);
		SeparatorTrack separator2 = new SeparatorTrack(Color.gray, 1, 0, Long.MAX_VALUE);
	
		CnaFlagTrack flag = new CnaFlagTrack(GBrowserConstants.BED_COLOR, Color.RED, 0, Long.MAX_VALUE);
		ScatterplotTrack freq = new ScatterplotTrack(GBrowserConstants.BED_COLOR, 0, Long.MAX_VALUE);
		ScatterplotTrack logRatio = new ScatterplotTrack(GBrowserConstants.BED_COLOR, 0, Long.MAX_VALUE);
				
		GBrowserView view = plot.getDataView();
		separator1.setView(view);
		separator2.setView(view);
		flag.setView(view);
		freq.setView(view);
		logRatio.setView(view);
		
		flag.setAreaRequestHandler(conversion);
		freq.setAreaRequestHandler(conversion);
		logRatio.setAreaRequestHandler(conversion);
		
		TrackGroup group = new TrackGroup(view);
		group.addTrack(flag);
		group.addTrack(separator1);
		group.addTrack(freq);
		group.addTrack(separator2);
		group.addTrack(logRatio);
		
		return group;
	}
}
