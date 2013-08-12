package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;

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
    protected CoverageTrack coverageTrack;
    protected CoverageAverageTrack coverageAverageTrack;
//    protected QualityCoverageTrack qualityCoverageTrack;
    protected DensityGraphTrack densityGraph;
    protected SeparatorTrack separatorReadPile;
    protected SeparatorTrack separatorReferenceSequence;
    protected SeparatorTrack separatorDensity;
    protected SeparatorTrack separatorQualityCoverage;
    protected SeparatorTrack sepTrackGel;
    protected SeparatorTrack separatorCoverageEstimate;
    
    private DataThread referenceSequenceFile;
	private DataThread details;
	private DataThread coverage;
	private DataThread estimate;
	
	private String title;
	private boolean initialised = false;


    public ReadTrackGroup(GBrowserView view, DataThread details, DataThread coverage, DataThread estimate,
    		DataThread seqFile, String title) {
        super(view);
        
        this.details = details;
        this.coverage = coverage;
        this.estimate = estimate;
        this.referenceSequenceFile = seqFile;
        this.title = title;    
    }

    public void initialise() {
        
        // Title
        titleTrack = new StatusTitleTrack(title, Color.black);  
        titleTrack.setView(view);
        tracks.add(titleTrack);
        
        if (details != null) {

        	// Detailed
        	readPileForward = new ReadPileTrack(referenceSequenceFile, fontColor);
        	readPileForward.setViewLimits(0, GBrowserConstants.SWITCH_VIEWS_AT);
        	readPileForward.setView(view);
        	readPileForward.addDataThread(details);
        	tracks.add(readPileForward);
        	separatorReadPile = new SeparatorTrack(Color.gray, 1);
        	separatorReadPile.setViewLimits(0, GBrowserConstants.SWITCH_VIEWS_AT);
        	separatorReadPile.setView(view);
        	separatorReadPile.setName("Reads");
        	tracks.add(separatorReadPile);
        	
        	titleTrack.addDataThread(details);
        }
        
        if (referenceSequenceFile != null) {
        	
            // Reference sequence		
            referenceSequence = new ReferenceSequenceTrack();
            referenceSequence.setViewLimits(0, GBrowserConstants.SHOW_REFERENCE_AT);
            referenceSequence.setView(view);
            referenceSequence.addDataThread(referenceSequenceFile);                        
            tracks.add(referenceSequence);            
            separatorReferenceSequence = new SeparatorTrack(Color.gray, 1);
            separatorReferenceSequence.setViewLimits(0, GBrowserConstants.SHOW_REFERENCE_AT);
            separatorReferenceSequence.setView(view);
            separatorReferenceSequence.setName("Reads");
            tracks.add(separatorReferenceSequence);
            
            titleTrack.addDataThread(referenceSequenceFile);
        }
        
        if (details != null) {
        	
        	// Overview
        	addCoverageEstimate();
	        
	        // Detailed - reversed
	        readsReverse = new ReadPileTrack(referenceSequenceFile, fontColor);
	        readsReverse.setViewLimits(0, GBrowserConstants.SWITCH_VIEWS_AT);
	        readsReverse.setView(view);
	        readsReverse.addDataThread(details);
	        readsReverse.setStrand(Strand.REVERSE);
	        tracks.add(readsReverse);
	    	SeparatorTrack sepTrackReads2 = new SeparatorTrack(Color.gray, 1);
	    	sepTrackReads2.setViewLimits(0, GBrowserConstants.SWITCH_VIEWS_AT);
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
	        
	        coverageTrack = new CoverageTrack(coverage, referenceSequenceFile);
	        coverageTrack.setViewLimits(0, GBrowserConstants.SHOW_AVERAGES);
	        coverageTrack.setView(view);
	        coverageTrack.setName("Coverage");
	        tracks.add(coverageTrack);
	        
	        coverageAverageTrack = new CoverageAverageTrack();
	        coverageAverageTrack.setViewLimits(GBrowserConstants.SHOW_AVERAGES, GBrowserConstants.SWITCH_VIEWS_AT);
	        coverageAverageTrack.addDataThread(coverage);
	        coverageAverageTrack.setView(view);
	        coverageAverageTrack.setName("Coverage");
	        tracks.add(coverageAverageTrack);
	        
	        titleTrack.addDataThread(coverage);
	        
	        separatorDensity = new SeparatorTrack(Color.gray, 1);
	        separatorDensity.setViewLimits(0, GBrowserConstants.SWITCH_VIEWS_AT);
	        separatorDensity.setView(view);
	        separatorDensity.setName("DensityGraphTrack");
	    	tracks.add(separatorDensity);	    
	
        	// Gel
        	densityGraph = new DensityGraphTrack(Color.WHITE);
        	densityGraph.setViewLimits(0, GBrowserConstants.SWITCH_VIEWS_AT);
        	densityGraph.setView(view);
        	densityGraph.addDataThread(coverage);
        	tracks.add(densityGraph);
        }
        
        this.initialised  = true;
    }

	protected void addCoverageEstimate() {
		coverageEstimate = new CoverageEstimateTrack();
		coverageEstimate.setViewLimits(GBrowserConstants.SWITCH_VIEWS_AT, Long.MAX_VALUE);
		coverageEstimate.setView(view);
		coverageEstimate.addDataThread(estimate);
		coverageEstimate.setName("CoverageEstimate");
		tracks.add(coverageEstimate);
		
		titleTrack.addDataThread(estimate);
	}
    
    public void setVisibleSNP(boolean b) {
    	check();
    	if (b) {
            readPileForward.enableSNPHighlight();
            readsReverse.enableSNPHighlight();
            coverageTrack.enableSNPHighlight();
        } else {
            readPileForward.disableSNPHiglight();
            readsReverse.disableSNPHiglight();
            coverageTrack.disableSNPHighlight();
        }
        view.fireDataRequests();
        view.redraw();
    }
    
    public void setStrandSpecificCoverageType(boolean b) {
    	
    	check();
    	
        coverageEstimate.setStrandSpecificCoverageType(b);
        coverageAverageTrack.setStrandSpecificCoverageType(b);
        coverageTrack.setStrandSpecificCoverageType(b);

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
