package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;

/**
 * Tracks containing information about reads: sequences themselves, gel,
 * profile etc.
 * 
 * @author Rimvydas Naktinis, Vilius Zukauskas
 *
 */
public class ReadTrackGroup extends TrackGroup {

	// Colors
//    private final Color histogramColor = new Color(100, 100, 100, 100); // translucent color
	private final Color histogramColor = new Color(100, 100, 100);
    private final Color fontColor = Color.black;

    // Tracks
    protected TitleTrack titleTrack;
    protected CoverageEstimateTrack readOverview;
    protected ReadPileTrack reads;
    protected ReferenceSequenceTrack seq;
    protected CoverageEstimateTrack readOverviewReversed;
    protected ReadPileTrack readsReversed;
    protected CoverageTrack profileTrack;
    protected CoverageTrack profileSNPTrack;
//    protected QualityCoverageTrack qualityCoverageTrack;
    protected GelTrack gelTrack;
    protected SeparatorTrack sepTrackReads;
    protected SeparatorTrack sepTrackSeq;
    protected SeparatorTrack sepTrackProfile;
    protected SeparatorTrack sepTrackProfileSNP;
    protected SeparatorTrack sepTrackQualityCoverage;
    protected SeparatorTrack sepTrackGel;
    
    private AreaRequestHandler seqFile;
	private AreaRequestHandler userData;
	private ReadpartDataProvider readpartProvider;
	private String title;
	private boolean initialised = false;
	private SeparatorTrack sepTrackReadOverview;

    public ReadTrackGroup(GBrowserView view, AreaRequestHandler userData,
    		AreaRequestHandler seqFile, String title) {
        super(view);
        
        this.userData = userData;
        this.seqFile = seqFile;
        this.title = title;
        
        if (userData != null) {
        	this.readpartProvider = new ReadpartDataProvider(view, userData);
        }
    }

    public void initialise() {
        
        // Title
        titleTrack = new TitleTrack(title, Color.black);
        tracks.add(titleTrack);
        
        if (userData != null) {
        	// Overview
        	addReadOverviewTrack();
        	sepTrackReadOverview = new SeparatorTrack(Color.gray, 1, GBrowserConstants.SWITCH_VIEWS_AT, Long.MAX_VALUE);
        	sepTrackReadOverview.setView(view);
        	tracks.add(sepTrackReadOverview);

        	// Detailed
        	reads = new ReadPileTrack(readpartProvider, seqFile, fontColor, 0, GBrowserConstants.SWITCH_VIEWS_AT);
        	reads.setView(view);
        	reads.setAreaRequestHandler(userData);
        	tracks.add(reads);
        	sepTrackReads = new SeparatorTrack(Color.gray, 1, 0, GBrowserConstants.SWITCH_VIEWS_AT);
        	sepTrackReads.setView(view);
        	sepTrackReads.setName("Reads");
        	tracks.add(sepTrackReads);
        }
        
        // Reference
        if (seqFile != null) {
            // Reference sequence		
            seq = new ReferenceSequenceTrack(GBrowserConstants.SHOW_REFERENCE_AT);
            seq.setView(view);
            seq.setAreaRequestHandler(seqFile);
                        
            tracks.add(seq);
            
            sepTrackSeq = new SeparatorTrack(Color.gray, 1, 0, GBrowserConstants.SHOW_REFERENCE_AT);
            sepTrackSeq.setView(view);
            sepTrackSeq.setName("Reads");
            tracks.add(sepTrackSeq);

        }
        
        if (userData != null) {
	        // Overview - reversed
	        addReadOverviewReversedTrack();
	        
	        // Detailed - reversed
	        readsReversed = new ReadPileTrack(readpartProvider, seqFile, fontColor, 0, GBrowserConstants.SWITCH_VIEWS_AT);
	        readsReversed.setView(view);
	        readsReversed.setAreaRequestHandler(userData);
	        readsReversed.setStrand(Strand.REVERSE);
	        tracks.add(readsReversed);
	    	SeparatorTrack sepTrackReads2 = new SeparatorTrack(Color.gray, 1, 0, GBrowserConstants.SWITCH_VIEWS_AT);
	    	sepTrackReads2.setView(view);
	    	sepTrackReads2.setName("Reads");
	        tracks.add(sepTrackReads2);
	        
	        // Profile
	        
	        Color forwardColor = new Color(0,0,0, 64);
	        Color reverseColor = new Color(
	        		GBrowserConstants.COLOR_BLUE.getRed(), 
	        		GBrowserConstants.COLOR_BLUE.getGreen(), 
	        		GBrowserConstants.COLOR_BLUE.getBlue(), 
	        		64);
	        Color totalColor = Color.gray;
	        
	//        profileTrack = new CoverageTrack(view, userData, readpartProvider, userDataHandler,
	//        		forwardColor, reverseColor, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT);
	//        profileTrack.setStrand(Strand.BOTH);
	//        tracks.add(profileTrack);
	//    	sepTrackProfile = new SeparatorTrack(view, Color.gray, 1, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT); 
	//    	sepTrackProfile.setName("ProfileTrack");
	//        tracks.add(sepTrackProfile);
	        
	        profileTrack = new CoverageTrack(readpartProvider, seqFile, forwardColor, reverseColor, 0, 
	        		GBrowserConstants.SWITCH_VIEWS_AT);
	        profileTrack.setView(view);
	        profileTrack.setAreaRequestHandler(userData);
	        profileTrack.setName("ProfileTrack");
	        tracks.add(profileTrack);
	        sepTrackProfile = new SeparatorTrack(Color.gray, 1, 0, GBrowserConstants.SWITCH_VIEWS_AT);
	        sepTrackProfile.setView(view);
	        sepTrackProfile.setName("ProfileTrack");
	    	tracks.add(sepTrackProfile);
	        
	        // SNP profile
	        profileSNPTrack = new CoverageTrack(readpartProvider, seqFile, totalColor, null, 0, 
	        		GBrowserConstants.SWITCH_VIEWS_AT);
	        profileSNPTrack.setView(view);
	        profileSNPTrack.setAreaRequestHandler(userData);
	        profileSNPTrack.setName("ProfileSNPTrack");
	        tracks.add(profileSNPTrack);
	    	sepTrackProfileSNP = new SeparatorTrack(Color.gray, 1, 0, GBrowserConstants.SWITCH_VIEWS_AT);
	    	sepTrackProfileSNP.setView(view);
	    	sepTrackProfileSNP.setName("ProfileSNPTrack");
	    	tracks.add(sepTrackProfileSNP);
	
        	// Gel
        	gelTrack = new GelTrack(readpartProvider, Color.WHITE, 0, GBrowserConstants.SWITCH_VIEWS_AT);
        	gelTrack.setView(view);
        	gelTrack.setAreaRequestHandler(userData);
        	gelTrack.setStrand(Strand.BOTH);
        	tracks.add(gelTrack);
        	//    	sepTrackGel = new SeparatorTrack(view, Color.gray, 1, 0, GBrowserConstants.SWITCH_VIEWS_AT); 
        	//    	sepTrackGel.setName("GelTrack");
        	//      tracks.add(sepTrackGel);
        }
        
        this.initialised  = true;
    }


	protected void addReadOverviewReversedTrack() {
		readOverviewReversed = new CoverageEstimateTrack(histogramColor, GBrowserConstants.SWITCH_VIEWS_AT, 
				false, true);
		readOverviewReversed.setView(view);
		readOverviewReversed.setAreaRequestHandler(userData);
        readOverviewReversed.setStrand(Strand.REVERSE);
        readOverviewReversed.setName("ReadOverview");
		tracks.add(readOverviewReversed);
	}

	protected void addReadOverviewTrack() {
		readOverview = new CoverageEstimateTrack(histogramColor, GBrowserConstants.SWITCH_VIEWS_AT, 
				false, true);
		readOverview.setView(view);
		readOverview.setAreaRequestHandler(userData);
		readOverview.setName("ReadOverview");
		tracks.add(readOverview);
	}
    
    public void setVisibleSNP(boolean b) {
    	check();
    	if (b) {
            reads.enableSNPHighlight(seqFile);
            readsReversed.enableSNPHighlight(seqFile);
            profileSNPTrack.enableSNPHighlight();
        } else {
            reads.disableSNPHiglight();
            readsReversed.disableSNPHiglight();
            profileSNPTrack.disableSNPHighlight();
        }
        view.fireAreaRequests();
        view.redraw();
    }


    @Override
    public String getName() {
    	check();
    	return "Read Track Group";
    }
    
    @Override
    public void showOrHide(String name, boolean state) {
    	check();
    	super.showOrHide(name, state);
    	if (name.equals("highlightSNP")) {
    		setVisibleSNP(state);
    	}
    }

	private void check() {
		if (!initialised) {
    		throw new IllegalStateException("you must call initialise() after creating this object");
    	}
	}
}
