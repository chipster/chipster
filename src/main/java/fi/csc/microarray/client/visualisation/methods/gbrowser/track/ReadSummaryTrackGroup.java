package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GenomeBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.TabixHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.constants.VisualConstants;

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
    protected CoverageTrack profileTrack;
    protected GelTrack gelTrack;
    protected SeqTrack seq;
    protected SeqBlockTrack readsReversed;
	private IntensityTrack readOverviewOld;
	private boolean hasReference = false;
	private SeparatorTrack sepTrackProfile;
	private Track profileSNPTrack;
	private SeparatorTrack sepTrackProfileSNP;
	private Track qualityCoverageTrack;
	private SeparatorTrack sepTrackQualityCoverage;
	private SeparatorTrack sepTrackGel;

    public ReadSummaryTrackGroup(View view, DataSource userData,
            Class<? extends AreaRequestHandler> userDataHandler, DataSource seqFile, File file) throws FileNotFoundException {
        super(view);
        
        Color histogramColor = Color.gray;
        Color fontColor = Color.black;
        
        // Title
        titleTrack = new TitleTrack(view, file.getName(), Color.black);

        // Overview - old
        readOverviewOld = new IntensityTrack(view, userData,
                userDataHandler, histogramColor, GenomeBrowserConstants.SWITCH_VIEWS_AT, false, true);

        // Overview
        readOverview = new TabixIntensityTrack(view, new TabixDataSource(file),
                TabixHandlerThread.class, histogramColor, SWITCH_VIEWS_AT, Long.MAX_VALUE);
        
        // Detailed
        reads = new SeqBlockTrack(view, userData,
                userDataHandler, fontColor, 0, SWITCH_VIEWS_AT);
        
        // Reference
        if (seqFile != null) {
        	hasReference = true;
            // Reference sequence
            seq = new SeqTrack(view, seqFile,
                    ChunkTreeHandlerThread.class, SHOW_REFERENCE_AT);
        }
        
        // Detailed
        readsReversed = new SeqBlockTrack(view, userData,
                userDataHandler, fontColor, 0, SWITCH_VIEWS_AT);
        readsReversed.setStrand(Strand.REVERSED);
        
        // Profile
        profileTrack = new CoverageTrack(view, userData, userDataHandler,
                Color.BLACK, VisualConstants.COLOR_BLUE, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT);
        profileTrack.setStrand(Strand.BOTH);
        
        // SNP profile
        profileSNPTrack = new CoverageAndSNPTrack(view, userData, userDataHandler, seqFile, ChunkTreeHandlerThread.class, 
                Color.BLACK, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT);
        profileSNPTrack.setStrand(Strand.BOTH); //Will be set anyway in the track constructor
        
        qualityCoverageTrack = new QualityCoverageTrack(view, userData, userDataHandler,
        		Color.ORANGE, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT);
        profileSNPTrack.setStrand(Strand.BOTH); //Will be set anyway in the track constructor
        
        // Gel
        gelTrack = new GelTrack(view, userData, userDataHandler,
                Color.WHITE, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT);
        gelTrack.setStrand(Strand.BOTH);
     
        // Add tracks to this group
        addTracks();
        
        this.setMenuVisible(false);
    }
    
    private void addTracks() {
        // Construct the list according to visibility
        this.tracks = new LinkedList<Track>();
        // Top separator
        tracks.add(titleTrack);
        tracks.add(readOverviewOld);
    	tracks.add(new SeparatorTrack(view, Color.gray, 1, (long)GenomeBrowserConstants.SWITCH_VIEWS_AT, Long.MAX_VALUE)); 
        tracks.add(readOverview);
        tracks.add(reads);
        tracks.add(new SeparatorTrack(view, Color.gray, 1, 0, Long.MAX_VALUE));
        
        // Only draw reference sequence if data is present
        if (hasReference ) {
            tracks.add(seq);
            tracks.add(new SeparatorTrack(view, Color.gray, 1, 0, SHOW_REFERENCE_AT));
        }

        tracks.add(readsReversed);

    	SeparatorTrack sepTrackReads2 = new SeparatorTrack(view, Color.gray, 1, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT); 
    	sepTrackReads2.setName("Reads");
        tracks.add(sepTrackReads2);

        tracks.add(profileTrack);
    	sepTrackProfile = new SeparatorTrack(view, Color.gray, 1, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT); 
    	sepTrackProfile.setName("ProfileTrack");
        tracks.add(sepTrackProfile);
        
        tracks.add(profileSNPTrack);
    	sepTrackProfileSNP = new SeparatorTrack(view, Color.gray, 1, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT);
    	sepTrackProfileSNP.setName("ProfileSNPTrack");
    	tracks.add(sepTrackProfileSNP);

        tracks.add(qualityCoverageTrack);
    	sepTrackQualityCoverage = new SeparatorTrack(view, Color.gray, 1, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT);
    	sepTrackQualityCoverage.setName("QualityCoverageTrack");
    	tracks.add(sepTrackQualityCoverage);
        
        tracks.add(gelTrack);
    	sepTrackGel = new SeparatorTrack(view, Color.gray, 1, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT); 
    	sepTrackGel.setName("GelTrack");
        tracks.add(sepTrackGel);
    }

    public void actionPerformed(ActionEvent e) {

    }
}
