package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.FileNotFoundException;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.view.View;

/**
 * Tracks containing information about reads: sequences themselves, gel,
 * profile etc.
 * 
 * @author Rimvydas Naktinis, Petri Klemelä
 *
 */
public class ReadSummaryTrackGroup extends ReadTrackGroup implements ActionListener {
        
    protected TabixIntensityTrack readOverviewSummary;
	private TabixDataSource summaryDataSource;

    public ReadSummaryTrackGroup(View view, DataSource userData, DataSource seqFile, String title, 
    		TabixDataSource summaryDataSource) throws FileNotFoundException {

    	super(view, userData, seqFile, title);
    	this.summaryDataSource = summaryDataSource;
    }

    @Override
    protected void addReadOverviewTrack() {
        readOverviewSummary = new TabixIntensityTrack(view, summaryDataSource, Color.black, 
        		GBrowserConstants.SWITCH_VIEWS_AT, Long.MAX_VALUE);
        tracks.add(readOverviewSummary);
    }

    @Override
    protected void addReadOverviewReversedTrack() {
    	// TODO add reversed summary overview track
    }

    
    public void actionPerformed(ActionEvent e) {

    }
}
