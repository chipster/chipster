package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JCheckBox;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TranscriptTrack.PartColor;

public class ReadTrackGroup extends TrackGroup implements ActionListener {
    
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
    protected GelTrack gelTrack;
    
    // Track switches
    private JCheckBox showGel = new JCheckBox("Gel track", true);
    private JCheckBox showProfile = new JCheckBox("Profile track", true);
    private JCheckBox showSNP = new JCheckBox("Highlight SNP", false);
    
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
        
        // Gel
        gelTrack = new GelTrack(view, userData, userDataHandler,
                Color.WHITE, 0, SWITCH_VIEWS_AT);
        gelTrack.setStrand(Strand.BOTH);
        
        // Add switches
        this.menu.addItem(showGel);
        this.menu.addItem(showProfile);
        this.menu.addItem(showSNP);
        //int startDrawing = 25;
        //showGel.setBounds(5, startDrawing, 100, 20);
        //showProfile.setBounds(5, startDrawing + 25, 100, 20);
        //showSNP.setBounds(5, startDrawing + 50, 100, 20);
        showGel.addActionListener(this);
        showProfile.addActionListener(this);
        showSNP.addActionListener(this);
        this.setMenuVisible(true);
    }
    
    @Override
    public List<Track> getTracks() {
        // Construct the list according to visibility
        List<Track> tracks = new LinkedList<Track>();
        tracks.add(titleTrack);
        tracks.add(readOverview);
        tracks.add(reads);
        tracks.add(new SeparatorTrack(view, Color.gray, 1, 0, Long.MAX_VALUE));
        
        // Only draw reference sequence if data is present
        if (hasReference) {
            tracks.add(seq);
            tracks.add(new SeparatorTrack(view, Color.gray, 1, 0, SHOW_REFERENCE_AT));
        }

        tracks.add(readOverviewReversed);
        tracks.add(readsReversed);
        
        // Only draw separator if profile track is visible
        if (showProfile.isSelected()) {
            tracks.add(new SeparatorTrack(view, Color.gray, 1, 0, SWITCH_VIEWS_AT));
            tracks.add(profileTrack);
        }
        
        // Only draw separator if gel track is visible
        if (showGel.isSelected()) {
            tracks.add(new SeparatorTrack(view, Color.gray, 1, 0, SWITCH_VIEWS_AT));
            tracks.add(gelTrack);
        }

        return tracks;
    }

    public void actionPerformed(ActionEvent e) {
        if (e.getSource() == showGel) {
            gelTrack.setVisible(showGel.isSelected());
            view.redraw();
        } else if (e.getSource() == showProfile) {
            profileTrack.setVisible(showProfile.isSelected());
            view.redraw();
        } else if (e.getSource() == showSNP && hasReference) {
            if (showSNP.isSelected()) {
                reads.enableSNPHighlight(seqFile, ChunkTreeHandlerThread.class);
                readsReversed.enableSNPHighlight(seqFile, ChunkTreeHandlerThread.class);
            } else {
                reads.disableSNPHiglight(seqFile);
                readsReversed.disableSNPHiglight(seqFile);
            }
            view.fireAreaRequests();
            view.redraw();
        }
    }

}
