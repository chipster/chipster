package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.TrackFactory;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.TabixHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TranscriptTrack.PartColor;

/**
 * Tracks containing information about reads: sequences themselves, gel,
 * profile etc.
 * 
 * @author naktinis, klemela
 *
 */
// FIXME integrate to ReadTrackGroup
public class ReadSummaryTrackGroup extends TrackGroup implements ActionListener {
        
    int SWITCH_VIEWS_AT = 50000;
    int SHOW_REFERENCE_AT = 800;

    // Tracks
    protected TitleTrack titleTrack;
    protected TabixIntensityTrack readOverview;
    protected SeqBlockTrack reads;
    protected ProfileTrack profileTrack;
    protected GelTrack gelTrack;
    protected SeqTrack seq;
//    protected IntensityTrack readOverviewReversed;
    protected SeqBlockTrack readsReversed;
    protected AcidProfile acidTrack;

    public ReadSummaryTrackGroup(View view, DataSource userData,
            Class<? extends AreaRequestHandler> userDataHandler, DataSource seqFile, String filename) throws FileNotFoundException {
        super(view);
        
        Color histogramColor = Color.gray;
        Color fontColor = Color.black;
        
        // Title
        titleTrack = new TitleTrack(view, filename, Color.black);
        
        // Overview
        readOverview = new TabixIntensityTrack(view, new TabixDataSource(new File(filename)),
                TabixHandlerThread.class, histogramColor, SWITCH_VIEWS_AT, Long.MAX_VALUE);
        
        // Detailed
        reads = new SeqBlockTrack(view, userData,
                userDataHandler, fontColor, 0, SWITCH_VIEWS_AT);
        
        // Reference
        if (seqFile != null) {
            // Reference sequence
//            hasReference = true;
//            this.seqFile = seqFile;
            seq = new SeqTrack(view, seqFile,
                    ChunkTreeHandlerThread.class, SHOW_REFERENCE_AT);
        }
        
//        // Overview
//        readOverviewReversed = new IntensityTrack(view, userData,
//                userDataHandler, histogramColor, SWITCH_VIEWS_AT);
//        readOverviewReversed.setStrand(Strand.REVERSED);
        
        // Detailed
        readsReversed = new SeqBlockTrack(view, userData,
                userDataHandler, fontColor, 0, SWITCH_VIEWS_AT);
        readsReversed.setStrand(Strand.REVERSED);
        
        // Profile
        profileTrack = new ProfileTrack(view, userData, userDataHandler,
                Color.BLACK, PartColor.CDS.c, 0, SWITCH_VIEWS_AT);
        profileTrack.setStrand(Strand.BOTH);
        
        // Acid profile
        acidTrack = new AcidProfile(view, userData, userDataHandler,
                0, SHOW_REFERENCE_AT);
        acidTrack.setStrand(Strand.FORWARD);        
        
        
        // Gel
        gelTrack = new GelTrack(view, userData, userDataHandler,
                Color.WHITE, 0, SWITCH_VIEWS_AT);
        gelTrack.setStrand(Strand.BOTH);
        
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
        tracks.add(reads);
        tracks.add(new SeparatorTrack(view, Color.gray, 1, 0, Long.MAX_VALUE));
        
        // Only draw reference sequence if data is present
//        if (hasReference) {
            tracks.add(seq);
            tracks.add(new SeparatorTrack(view, Color.gray, 1, 0, SHOW_REFERENCE_AT));
//        }

//        tracks.add(readOverviewReversed);
        tracks.add(readsReversed);
        
        // Only draw separator if profile track is visible
//        if (showProfile.isSelected()) {
            tracks.add(new SeparatorTrack(view, Color.gray, 1, 0, SWITCH_VIEWS_AT));
            tracks.add(profileTrack);
//        }
        
//        if (showAcid.isSelected()) {
        	tracks.add(new SeparatorTrack(view, Color.gray, 1, 0, SWITCH_VIEWS_AT));
            tracks.add(acidTrack);
//        
        
        // Only draw separator if gel track is visible
//        if (showGel.isSelected()) {
            tracks.add(new SeparatorTrack(view, Color.gray, 1, 0, SWITCH_VIEWS_AT));
            tracks.add(gelTrack);
//        }
    }

    public void actionPerformed(ActionEvent e) {

    }
}
