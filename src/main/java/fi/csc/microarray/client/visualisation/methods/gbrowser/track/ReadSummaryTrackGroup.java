package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.FileNotFoundException;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;

/**
 * Tracks containing information about reads: sequences themselves, gel,
 * profile etc.
 * 
 * @author Rimvydas Naktinis, Petri Klemelä
 *
 */
public class ReadSummaryTrackGroup extends ReadTrackGroup implements ActionListener {
        
    protected TabixIntensityTrack readOverviewSummary;
	private AreaRequestHandler summaryDataSource;

    public ReadSummaryTrackGroup(GBrowserView view, AreaRequestHandler details, AreaRequestHandler estimate, AreaRequestHandler seqFile, String title, 
    		AreaRequestHandler summaryDataSource) throws FileNotFoundException {

    	super(view, details, estimate, seqFile, title);
    	this.summaryDataSource = summaryDataSource;
    }

    @Override
    protected void addCoverageEstimate() {
        readOverviewSummary = new TabixIntensityTrack(Color.black, 
        		GBrowserConstants.SWITCH_VIEWS_AT, Long.MAX_VALUE);
        readOverviewSummary.setView(view);
        readOverviewSummary.addAreaRequestHandler(summaryDataSource);
        tracks.add(readOverviewSummary);
    }
    
    public void actionPerformed(ActionEvent e) {

    }
}
