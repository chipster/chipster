package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TranscriptTrack.PartColor;

public class ReadTrackGroup extends TrackGroup {
    
    // Constants
    int SWITCH_VIEWS_AT = 50000;
    int SHOW_REFERENCE_AT = 800;
    
    // Tracks
    protected IntensityTrack readOverview;
    protected SeqBlockTrack reads;
    protected SeqTrack seq;
    protected IntensityTrack readOverviewReversed;
    protected SeqBlockTrack readsReversed;
    protected ProfileTrack profileTrack;
    protected GelTrack gelTrack;
    
    private boolean hasReference = false;

    public ReadTrackGroup(View view, DataSource userData,
            Class<? extends AreaRequestHandler> userDataHandler,
            DataSource seqFile) {
        super(view);
        
        Color histogramColor = Color.gray;
        Color fontColor = Color.black;
        
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
        
        // Gel
        gelTrack = new GelTrack(view, userData, userDataHandler,
                Color.WHITE, 0, SWITCH_VIEWS_AT);
        gelTrack.setStrand(Strand.BOTH);
    }
    
    @Override
    public List<Track> getTracks() {
        // Construct the list according to visibility
        List<Track> tracks = new LinkedList<Track>();
        tracks.add(readOverview);
        tracks.add(reads);
        tracks.add(new SeparatorTrack(view, Color.gray, 1, 0, Long.MAX_VALUE));
        
        if (hasReference) {
            tracks.add(seq);
            tracks.add(new SeparatorTrack(view, Color.gray, 1, 0, SHOW_REFERENCE_AT));
        }

        tracks.add(readOverviewReversed);
        tracks.add(readsReversed);
        tracks.add(new SeparatorTrack(view, Color.gray, 1, 0, SWITCH_VIEWS_AT));
        
        tracks.add(profileTrack);
        tracks.add(new SeparatorTrack(view, Color.gray, 1, 0, SWITCH_VIEWS_AT));
        tracks.add(gelTrack);

        return tracks;
    }

}
