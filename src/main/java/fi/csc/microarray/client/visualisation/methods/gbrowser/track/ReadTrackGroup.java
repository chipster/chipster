package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
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
    protected IntensityTrack readOverview;
    protected ReadPileTrack reads;
    protected ReferenceSequenceTrack seq;
    protected IntensityTrack readOverviewReversed;
    protected ReadPileTrack readsReversed;
    protected CoverageAndSNPTrack profileTrack;
    protected CoverageAndSNPTrack profileSNPTrack;
//    protected QualityCoverageTrack qualityCoverageTrack;
    protected GelTrack gelTrack;
    protected SeparatorTrack sepTrackReads;
    protected SeparatorTrack sepTrackSeq;
    protected SeparatorTrack sepTrackProfile;
    protected SeparatorTrack sepTrackProfileSNP;
    protected SeparatorTrack sepTrackQualityCoverage;
    protected SeparatorTrack sepTrackGel;
    
    private DataSource seqFile;
	private DataSource userData;
	private ReadpartDataProvider readpartProvider;
	private String title;
	private boolean initialised = false;
	private SeparatorTrack sepTrackReadOverview;

    public ReadTrackGroup(GBrowserView view, DataSource userData,
            DataSource seqFile, String title) {
        super(view);
        
        this.userData = userData;
        this.readpartProvider = new ReadpartDataProvider(view, userData);
        this.seqFile = seqFile;
        this.title = title;
    }

    public void initialise() {
        
        // Title
        titleTrack = new TitleTrack(view, title, Color.black);
        tracks.add(titleTrack);
        
        // Overview
        addReadOverviewTrack();
        sepTrackReadOverview = new SeparatorTrack(view, Color.gray, 1, GBrowserConstants.SWITCH_VIEWS_AT, Long.MAX_VALUE);
        tracks.add(sepTrackReadOverview);

        // Detailed
        reads = new ReadPileTrack(view, userData, readpartProvider, fontColor, 0, GBrowserConstants.SWITCH_VIEWS_AT);
        tracks.add(reads);
        sepTrackReads = new SeparatorTrack(view, Color.gray, 1, 0, GBrowserConstants.SWITCH_VIEWS_AT);
        sepTrackReads.setName("Reads");
        tracks.add(sepTrackReads);
        
        // Reference
        if (seqFile != null) {
            // Reference sequence		
            seq = new ReferenceSequenceTrack(view, seqFile, GBrowserConstants.SHOW_REFERENCE_AT);
                        
            tracks.add(seq);
            
            sepTrackSeq = new SeparatorTrack(view, Color.gray, 1, 0, GBrowserConstants.SHOW_REFERENCE_AT);
            sepTrackSeq.setName("Reads");
            tracks.add(sepTrackSeq);

        }
        
        // Overview - reversed
        addReadOverviewReversedTrack();
        
        // Detailed - reversed
        readsReversed = new ReadPileTrack(view, userData, readpartProvider, fontColor, 0, GBrowserConstants.SWITCH_VIEWS_AT);
        readsReversed.setStrand(Strand.REVERSED);
        tracks.add(readsReversed);
    	SeparatorTrack sepTrackReads2 = new SeparatorTrack(view, Color.gray, 1, 0, GBrowserConstants.SWITCH_VIEWS_AT); 
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
        
        profileTrack = new CoverageAndSNPTrack(view, userData, readpartProvider, seqFile, forwardColor, reverseColor, 0, 
        		GBrowserConstants.SWITCH_VIEWS_AT);
        profileTrack.setName("ProfileTrack");
        tracks.add(profileTrack);
        sepTrackProfile = new SeparatorTrack(view, Color.gray, 1, 0, GBrowserConstants.SWITCH_VIEWS_AT);
        sepTrackProfile.setName("ProfileTrack");
    	tracks.add(sepTrackProfile);
        
        // SNP profile
        profileSNPTrack = new CoverageAndSNPTrack(view, userData, readpartProvider, seqFile, totalColor, null, 0, 
        		GBrowserConstants.SWITCH_VIEWS_AT);
        profileSNPTrack.setName("ProfileSNPTrack");
        tracks.add(profileSNPTrack);
    	sepTrackProfileSNP = new SeparatorTrack(view, Color.gray, 1, 0, GBrowserConstants.SWITCH_VIEWS_AT);
    	sepTrackProfileSNP.setName("ProfileSNPTrack");
    	tracks.add(sepTrackProfileSNP);

    	// Quality coverage
//        qualityCoverageTrack = new QualityCoverageTrack(view, userData, userDataHandler,
//        		Color.ORANGE, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT);
//        tracks.add(qualityCoverageTrack);
//    	sepTrackQualityCoverage = new SeparatorTrack(view, Color.gray, 1, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT);
//    	sepTrackQualityCoverage.setName("QualityCoverageTrack");
//    	tracks.add(sepTrackQualityCoverage);
        
        // Gel
        gelTrack = new GelTrack(view, userData, readpartProvider, Color.WHITE, 0, GBrowserConstants.SWITCH_VIEWS_AT);
        gelTrack.setStrand(Strand.BOTH);
        tracks.add(gelTrack);
    	sepTrackGel = new SeparatorTrack(view, Color.gray, 1, 0, GBrowserConstants.SWITCH_VIEWS_AT); 
    	sepTrackGel.setName("GelTrack");
        tracks.add(sepTrackGel);
        
        this.initialised  = true;
    }


	protected void addReadOverviewReversedTrack() {
		readOverviewReversed = new IntensityTrack(view, userData, histogramColor, GBrowserConstants.SWITCH_VIEWS_AT, 
				false, true);
        readOverviewReversed.setStrand(Strand.REVERSED);
        readOverviewReversed.setName("ReadOverview");
		tracks.add(readOverviewReversed);
	}

	protected void addReadOverviewTrack() {
		readOverview = new IntensityTrack(view, userData, histogramColor, GBrowserConstants.SWITCH_VIEWS_AT, 
				false, true);
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
	
	@Override
	public boolean isFixedHeight() {
		return false;
	}
}
