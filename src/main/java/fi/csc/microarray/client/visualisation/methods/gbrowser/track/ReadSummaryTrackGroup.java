package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.TrackFactory;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;

/**
 * Tracks containing information about reads: sequences themselves, gel,
 * profile etc.
 * 
 * @author naktinis, klemela
 *
 */
public class ReadSummaryTrackGroup extends TrackGroup implements ActionListener {
        
    // Tracks
    protected TitleTrack titleTrack;
    protected TabixIntensityTrack readOverview;
    protected SeqBlockTrack reads;
    protected ProfileTrack profileTrack;
    protected GelTrack gelTrack;
            
    public ReadSummaryTrackGroup(View view, DataSource userData,
            Class<? extends AreaRequestHandler> userDataHandler, String title) {
        super(view);
        
        Color histogramColor = Color.gray;
//        Color fontColor = Color.black;
        
        // Title
        titleTrack = new TitleTrack(view, title, Color.black);
        
        // Overview
        readOverview = new TabixIntensityTrack(view, userData,
                userDataHandler, histogramColor, 0, Long.MAX_VALUE);
        
//        // Detailed
//        reads = new SeqBlockTrack(view, userData,
//                userDataHandler, fontColor, 0, Long.MAX_VALUE);      
//        
//        // Profile
//        profileTrack = new ProfileTrack(view, userData, userDataHandler,
//                Color.BLACK, PartColor.CDS.c, 0, Long.MAX_VALUE);
//        profileTrack.setStrand(Strand.BOTH);
//        
//        // Gel
//        gelTrack = new GelTrack(view, userData, userDataHandler,
//                Color.WHITE, 0, Long.MAX_VALUE);
//        gelTrack.setStrand(Strand.BOTH);
        
        // Add tracks to this group
        addTracks();
        
        this.setMenuVisible(false);
    }
    
    private void addTracks() {
        // Construct the list according to visibility
        this.tracks = new LinkedList<Track>();
        // Top separator
        tracks.add(TrackFactory.createThickSeparatorTrack(view));
        tracks.add(titleTrack);
        tracks.add(readOverview);
//        tracks.add(reads);
//        tracks.add(new SeparatorTrack(view, Color.gray, 1, 0, Long.MAX_VALUE));
//        tracks.add(new SeparatorTrack(view, Color.gray, 1, 0, Long.MAX_VALUE));
//        tracks.add(profileTrack);
//        tracks.add(new SeparatorTrack(view, Color.gray, 1, 0, Long.MAX_VALUE));
//        tracks.add(gelTrack);
    }

    public void actionPerformed(ActionEvent e) {

    }
}
