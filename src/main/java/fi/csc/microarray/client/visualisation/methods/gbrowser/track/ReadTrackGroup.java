package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.LinkedList;

import javax.swing.JCheckBox;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.TrackFactory;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TranscriptTrack.PartColor;

/**
 * Tracks containing information about reads: sequences themselves, gel,
 * profile etc.
 * 
 * @author naktinis, zukauska
 *
 */
public class ReadTrackGroup extends TrackGroup {
    
    // Constants
    int SWITCH_VIEWS_AT = 50000;
    int SHOW_REFERENCE_AT = 800;
    
    // Tracks
    protected TitleTrack titleTrack;
    protected IntensityTrack readOverview;
    protected SeqBlockTrack reads;
    protected SeqTrack seq;
    protected IntensityTrack readOverviewReversed;
    protected SeqBlockTrack readsReversed;
    protected ProfileTrack profileTrack;
    protected AcidProfile acidTrack;
    protected GelTrack gelTrack;
    protected Track sepTrackTitle;
    protected SeparatorTrack sepTrackReads;
    protected SeparatorTrack sepTrackSeq;
    protected SeparatorTrack sepTrackProfile;
    protected SeparatorTrack sepTrackAcid;
    protected SeparatorTrack sepTrackGel;
    
    // Reference sequence
    private DataSource seqFile;
    private boolean hasReference = false;

    public ReadTrackGroup(View view, DataSource userData,
            Class<? extends AreaRequestHandler> userDataHandler,
            DataSource seqFile, String title) {
        super(view);
        
        Color histogramColor = Color.gray;
        Color fontColor = Color.black;
        
        // Title
        titleTrack = new TitleTrack(view, title, Color.black);
        
        // Overview
        readOverview = new IntensityTrack(view, userData,
                userDataHandler, histogramColor, SWITCH_VIEWS_AT);
        
        // Detailed
        reads = new SeqBlockTrack(view, userData,
                userDataHandler, fontColor, 0, SWITCH_VIEWS_AT);
        
        // Reference
        if (seqFile != null) {
            // Reference sequence
            hasReference = true;
            this.seqFile = seqFile;
            seq = new SeqTrack(view, seqFile,
                    ChunkTreeHandlerThread.class, SHOW_REFERENCE_AT);
        }
        
        // Overview
        readOverviewReversed = new IntensityTrack(view, userData,
                userDataHandler, histogramColor, SWITCH_VIEWS_AT);
        readOverviewReversed.setStrand(Strand.REVERSED);
        
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

    }
    
    private void addTracks() {
        // Construct the list according to visibility
        this.tracks = new LinkedList<Track>();
        // Top separator
        sepTrackTitle = TrackFactory.createThickSeparatorTrack(view); 
        tracks.add(sepTrackTitle);
        tracks.add(titleTrack);
        tracks.add(readOverview);
        tracks.add(reads);
        sepTrackReads = new SeparatorTrack(view, Color.gray, 1, 0, Long.MAX_VALUE); 
        tracks.add(sepTrackReads);
        
        // Only draw reference sequence if data is present
        if (hasReference) {
            tracks.add(seq);
            sepTrackSeq = new SeparatorTrack(view, Color.gray, 1, 0, SHOW_REFERENCE_AT); 
            tracks.add(sepTrackSeq);
        }

        tracks.add(readOverviewReversed);
        tracks.add(readsReversed);
        
        // Only draw separator if profile track is visible
    	sepTrackProfile = new SeparatorTrack(view, Color.gray, 1, 0, SWITCH_VIEWS_AT); 
        tracks.add(sepTrackProfile);
        tracks.add(profileTrack);
        
    	sepTrackAcid = new SeparatorTrack(view, Color.gray, 1, 0, SWITCH_VIEWS_AT); 
    	tracks.add(sepTrackAcid);
        tracks.add(acidTrack);
        
        // Only draw separator if gel track is visible
    	sepTrackGel = new SeparatorTrack(view, Color.gray, 1, 0, SWITCH_VIEWS_AT); 
        tracks.add(sepTrackGel);
        tracks.add(gelTrack);
    }
    
    public void setVisibleSNP(boolean b) {
    	if (b) {
            reads.enableSNPHighlight(seqFile, ChunkTreeHandlerThread.class);
            readsReversed.enableSNPHighlight(seqFile, ChunkTreeHandlerThread.class);
        } else {
            reads.disableSNPHiglight(seqFile);
            readsReversed.disableSNPHiglight(seqFile);
        }
        view.fireAreaRequests();
        view.redraw();
    }


    @Override
    public String getName() {
    	return "Read Track Group";
    }
    
    @Override
    public void showOrHide(String name, boolean state) {
    	super.showOrHide(name, state);
    	if (name.equals("highlightSNP")) {
    		setVisibleSNP(state);
    	}
    }
    
}
