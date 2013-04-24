package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;

/**
 * Tracks containing information about reads: sequences themselves, gel,
 * profile etc.
 * 
 * @author Rimvydas Naktinis, Vilius Zukauskas
 *
 */
public class ReadTrackGroup extends TrackGroup {

	// Colors
    private final Color fontColor = Color.black;

    // Tracks
    protected StatusTitleTrack titleTrack;
    protected CoverageEstimateTrack coverageEstimate;
    protected ReadPileTrack readPileForward;
    protected ReferenceSequenceTrack referenceSequence;
    
    protected ReadPileTrack readsReverse;
    protected CoverageTrack coverage;
//    protected QualityCoverageTrack qualityCoverageTrack;
    protected GelTrack gelTrack;
    protected SeparatorTrack separatorReadPile;
    protected SeparatorTrack separatorReferenceSequence;
    protected SeparatorTrack separatorCoverage;
    protected SeparatorTrack separatorQualityCoverage;
    protected SeparatorTrack sepTrackGel;
    protected SeparatorTrack separatorCoverageEstimate;
    
    private AreaRequestHandler referenceSequenceFile;
	private AreaRequestHandler details;
	private AreaRequestHandler estimate;
	private ReadpartDataProvider readpartProvider;
	private String title;
	private boolean initialised = false;

    public ReadTrackGroup(GBrowserView view, AreaRequestHandler details, AreaRequestHandler estimate,
    		AreaRequestHandler seqFile, String title) {
        super(view);
        
        this.details = details;
        this.estimate = estimate;
        this.referenceSequenceFile = seqFile;
        this.title = title;
        
        if (details != null) {
        	this.readpartProvider = new ReadpartDataProvider(view, details);
        }
    }

    public void initialise() {
        
        // Title
        titleTrack = new StatusTitleTrack(title, Color.black);  
        titleTrack.setView(view);
        tracks.add(titleTrack);
        
        if (details != null) {

        	// Detailed
        	readPileForward = new ReadPileTrack(readpartProvider, referenceSequenceFile, fontColor, 0, GBrowserConstants.SWITCH_VIEWS_AT);
        	readPileForward.setView(view);
        	readPileForward.addAreaRequestHandler(details);
        	tracks.add(readPileForward);
        	separatorReadPile = new SeparatorTrack(Color.gray, 1, 0, GBrowserConstants.SWITCH_VIEWS_AT);
        	separatorReadPile.setView(view);
        	separatorReadPile.setName("Reads");
        	tracks.add(separatorReadPile);
        	
        	titleTrack.addAreaRequestHandler(details);
        }
        
        if (referenceSequenceFile != null) {
        	
            // Reference sequence		
            referenceSequence = new ReferenceSequenceTrack(GBrowserConstants.SHOW_REFERENCE_AT);
            referenceSequence.setView(view);
            referenceSequence.addAreaRequestHandler(referenceSequenceFile);                        
            tracks.add(referenceSequence);            
            separatorReferenceSequence = new SeparatorTrack(Color.gray, 1, 0, GBrowserConstants.SHOW_REFERENCE_AT);
            separatorReferenceSequence.setView(view);
            separatorReferenceSequence.setName("Reads");
            tracks.add(separatorReferenceSequence);
            
            titleTrack.addAreaRequestHandler(referenceSequenceFile);
        }
        
        if (details != null) {
        	
        	// Overview
        	addCoverageEstimate();
	        
	        // Detailed - reversed
	        readsReverse = new ReadPileTrack(readpartProvider, referenceSequenceFile, fontColor, 0, GBrowserConstants.SWITCH_VIEWS_AT);
	        readsReverse.setView(view);
	        readsReverse.addAreaRequestHandler(details);
	        readsReverse.setStrand(Strand.REVERSE);
	        tracks.add(readsReverse);
	    	SeparatorTrack sepTrackReads2 = new SeparatorTrack(Color.gray, 1, 0, GBrowserConstants.SWITCH_VIEWS_AT);
	    	sepTrackReads2.setView(view);
	    	sepTrackReads2.setName("Reads");
	        tracks.add(sepTrackReads2);
	        
	        // Profile	   
	        
	//        profileTrack = new CoverageTrack(view, userData, readpartProvider, userDataHandler,
	//        		forwardColor, reverseColor, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT);
	//        profileTrack.setStrand(Strand.BOTH);
	//        tracks.add(profileTrack);
	//    	sepTrackProfile = new SeparatorTrack(view, Color.gray, 1, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT); 
	//    	sepTrackProfile.setName("ProfileTrack");
	//        tracks.add(sepTrackProfile);
	        
	        coverage = new CoverageTrack(readpartProvider, details, referenceSequenceFile, 0, GBrowserConstants.SWITCH_VIEWS_AT);
	        coverage.setView(view);
	        coverage.setName("Coverage");
	        tracks.add(coverage);
	        separatorCoverage = new SeparatorTrack(Color.gray, 1, 0, GBrowserConstants.SWITCH_VIEWS_AT);
	        separatorCoverage.setView(view);
	        separatorCoverage.setName("Coverage");
	    	tracks.add(separatorCoverage);	    
	
        	// Gel
        	gelTrack = new GelTrack(readpartProvider, Color.WHITE, 0, GBrowserConstants.SWITCH_VIEWS_AT);
        	gelTrack.setView(view);
        	gelTrack.addAreaRequestHandler(details);
        	gelTrack.setStrand(Strand.BOTH);
        	tracks.add(gelTrack);
        	//    	sepTrackGel = new SeparatorTrack(view, Color.gray, 1, 0, GBrowserConstants.SWITCH_VIEWS_AT); 
        	//    	sepTrackGel.setName("GelTrack");
        	//      tracks.add(sepTrackGel);
        }
        
        this.initialised  = true;
    }

	protected void addCoverageEstimate() {
		coverageEstimate = new CoverageEstimateTrack(GBrowserConstants.SWITCH_VIEWS_AT);
		coverageEstimate.setView(view);
		coverageEstimate.addAreaRequestHandler(estimate);
		coverageEstimate.setName("CoverageEstimate");
		tracks.add(coverageEstimate);
		
		titleTrack.addAreaRequestHandler(estimate);
	}
    
    public void setVisibleSNP(boolean b) {
    	check();
    	if (b) {
            readPileForward.enableSNPHighlight();
            readsReverse.enableSNPHighlight();
            coverage.enableSNPHighlight();
        } else {
            readPileForward.disableSNPHiglight();
            readsReverse.disableSNPHiglight();
            coverage.disableSNPHighlight();
        }
        view.fireAreaRequests();
        view.redraw();
    }
    
    public void setStrandSpecificCoverageType(boolean b) {
    	
    	check();
    	
        coverageEstimate.setStrandSpecificCoverageType(b);
        coverage.setStrandSpecificCoverageType(b);

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
    		
    	} else if (name.equals("StrandSpecificCoverageType")) {
    		setStrandSpecificCoverageType(state);
    	}
    }

	private void check() {
		if (!initialised) {
    		throw new IllegalStateException("you must call initialise() after creating this object");
    	}
	}
}
