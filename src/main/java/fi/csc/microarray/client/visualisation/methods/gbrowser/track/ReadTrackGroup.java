package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GenomeBrowserConstants;
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

	// Colors
//    private final Color histogramColor = new Color(100, 100, 100, 100); // translucent color
	private final Color histogramColor = new Color(100, 100, 100);
    private final Color fontColor = Color.black;

    // Tracks
    protected TitleTrack titleTrack;
    protected IntensityTrack readOverview;
    protected SeqBlockTrack reads;
    protected SeqTrack seq;
    protected IntensityTrack readOverviewReversed;
    protected SeqBlockTrack readsReversed;
    protected CoverageTrack profileTrack;
    protected CoverageAndSNPTrack profileSNPTrack;
    protected QualityCoverageTrack qualityCoverageTrack;
    protected GelTrack gelTrack;
    protected SeparatorTrack sepTrackReads;
    protected SeparatorTrack sepTrackSeq;
    protected SeparatorTrack sepTrackProfile;
    protected SeparatorTrack sepTrackProfileSNP;
    protected SeparatorTrack sepTrackQualityCoverage;
    protected SeparatorTrack sepTrackGel;
    
    // Reference sequence
    private DataSource seqFile;
    private boolean hasReference = false;

    public ReadTrackGroup(View view, DataSource userData,
            Class<? extends AreaRequestHandler> userDataHandler,
            DataSource seqFile, String title) {
        super(view);
        
        // Title
        titleTrack = new TitleTrack(view, title, Color.black);
        
        // Overview
        readOverview = new IntensityTrack(view, userData,
                userDataHandler, histogramColor, GenomeBrowserConstants.SWITCH_VIEWS_AT, false);
            
        reads = new SeqBlockTrack(view, userData,
                userDataHandler, fontColor, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT);
        
        // Reference
        if (seqFile != null) {
            // Reference sequence
            hasReference = true;
            this.seqFile = seqFile;
            seq = new SeqTrack(view, seqFile,
                    ChunkTreeHandlerThread.class, GenomeBrowserConstants.SHOW_REFERENCE_AT);
        }
        
        // Overview
        readOverviewReversed = new IntensityTrack(view, userData,
                userDataHandler, histogramColor, GenomeBrowserConstants.SWITCH_VIEWS_AT, false);
        readOverviewReversed.setStrand(Strand.REVERSED);
        
        // Detailed
        readsReversed = new SeqBlockTrack(view, userData,
                userDataHandler, fontColor, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT);
        readsReversed.setStrand(Strand.REVERSED);
        
        // Profile
        profileTrack = new CoverageTrack(view, userData, userDataHandler,
                Color.BLACK, PartColor.CDS.c, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT);
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
    }
    
    private void addTracks() {

    	this.tracks = new LinkedList<Track>();
        tracks.add(titleTrack);
        tracks.add(readOverview);
        tracks.add(reads);
        sepTrackReads = new SeparatorTrack(view, Color.gray, 1, 0, Long.MAX_VALUE);
        sepTrackReads.setName("Reads");
        tracks.add(sepTrackReads);
        
        // Only draw reference sequence if data is present
        if (hasReference) {
            tracks.add(seq);
            sepTrackSeq = new SeparatorTrack(view, Color.gray, 1, 0, GenomeBrowserConstants.SHOW_REFERENCE_AT);
            sepTrackSeq.setName("Reads");
            tracks.add(sepTrackSeq);
        }

        tracks.add(readOverviewReversed);
        tracks.add(readsReversed);
        
    	sepTrackProfile = new SeparatorTrack(view, Color.gray, 1, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT); 
    	sepTrackProfile.equals("ProfileTrack");
        tracks.add(sepTrackProfile);
        tracks.add(profileTrack);
        
    	sepTrackProfileSNP = new SeparatorTrack(view, Color.gray, 1, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT);
    	sepTrackProfileSNP.setName("ProfileSNPTrack");
    	tracks.add(sepTrackProfileSNP);
        tracks.add(profileSNPTrack);

    	sepTrackQualityCoverage = new SeparatorTrack(view, Color.gray, 1, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT);
    	sepTrackQualityCoverage.setName("QualityCoverageTrack");
    	tracks.add(sepTrackQualityCoverage);
        tracks.add(qualityCoverageTrack);
        
    	sepTrackGel = new SeparatorTrack(view, Color.gray, 1, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT); 
    	sepTrackGel.setName("GelTrack");
        tracks.add(sepTrackGel);
        tracks.add(gelTrack);
    }
    
    public void setVisibleSNP(boolean b) {
    	if (b) {
            reads.enableSNPHighlight(seqFile, ChunkTreeHandlerThread.class);
            readsReversed.enableSNPHighlight(seqFile, ChunkTreeHandlerThread.class);
            profileSNPTrack.enableSNPHighlight();
        } else {
            reads.disableSNPHiglight(seqFile);
            readsReversed.disableSNPHiglight(seqFile);
            profileSNPTrack.disableSNPHighlight();
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
