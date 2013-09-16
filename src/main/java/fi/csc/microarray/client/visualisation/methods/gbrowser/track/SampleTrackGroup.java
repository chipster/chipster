package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserSettings.CoverageType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LayoutTool.LayoutMode;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;

/**
 * Tracks containing information about reads: sequences themselves, gel,
 * profile etc.
 * 
 * @author Rimvydas Naktinis, Vilius Zukauskas, Petri Klemelä
 *
 */
public class SampleTrackGroup extends TrackGroup {

	// Colors
    private final Color fontColor = Color.black;

    // Tracks
    protected CoverageEstimateTrack coverageEstimateTrack;
    protected ReadPileTrack readPileForward;
    protected ReferenceSequenceTrack referenceSequence;
    
    protected ReadPileTrack readPileReverse;
    protected CoverageTrack coverageTrack;
    protected CoverageAverageTrack coverageAverageTrack;
//    protected QualityCoverageTrack qualityCoverageTrack;
    protected DensityGraphTrack densityGraphTrack;
    protected SeparatorTrack separatorReadPile;
    protected SeparatorTrack separatorReferenceSequence;
    protected SeparatorTrack separatorDensity;
    protected SeparatorTrack separatorQualityCoverage;
    protected SeparatorTrack sepTrackGel;
    protected SeparatorTrack separatorCoverageEstimate;
    
    private DataThread referenceSequenceFile;
	private DataThread detailsDataThread;
	private DataThread coverageDataThread;
	private DataThread estimateDataThread;

	private boolean strandSpecific;

	private boolean coverage;

	private boolean coverageEstimate;

	private boolean reads;

	private boolean highlightSnp;

	private boolean densityGraph;

	private boolean fullMode;


    public SampleTrackGroup(GBrowserView view, DataThread details, DataThread coverage, DataThread estimate,
    		DataThread seqFile, String title) {
        super(view);
        
        super.setName(title);
        
        this.detailsDataThread = details;
        this.coverageDataThread = coverage;
        this.estimateDataThread = estimate;
        this.referenceSequenceFile = seqFile;
        
        setSettingsEnabled(true);
    }

    public void initialise() {
        
        if (detailsDataThread != null) {

        	if (reads) {
        		// Detailed
        		readPileForward = new ReadPileTrack(referenceSequenceFile, fontColor);
        		readPileForward.setViewLimits(0, GBrowserConstants.SWITCH_VIEWS_AT);
        		readPileForward.setView(view);
        		readPileForward.addDataThread(detailsDataThread);
                readPileForward.setSNPHighlight(highlightSnp);
        		addTrack(readPileForward);
        		separatorReadPile = new SeparatorTrack(Color.gray, 1);
        		separatorReadPile.setViewLimits(0, GBrowserConstants.SWITCH_VIEWS_AT);
        		separatorReadPile.setView(view);
        		separatorReadPile.setTrackName("Reads");
        		addTrack(separatorReadPile);

        		getStatusAnimation().addDataThread(detailsDataThread);
        	}
        }
        
        if (referenceSequenceFile != null) {
        	
            // Reference sequence		
            referenceSequence = new ReferenceSequenceTrack();
            referenceSequence.setViewLimits(0, GBrowserConstants.SHOW_REFERENCE_AT);
            referenceSequence.setView(view);
            referenceSequence.addDataThread(referenceSequenceFile);                        
            addTrack(referenceSequence);            
            separatorReferenceSequence = new SeparatorTrack(Color.gray, 1);
            separatorReferenceSequence.setViewLimits(0, GBrowserConstants.SHOW_REFERENCE_AT);
            separatorReferenceSequence.setView(view);
            separatorReferenceSequence.setTrackName("Reads");
            addTrack(separatorReferenceSequence);
            
            getStatusAnimation().addDataThread(referenceSequenceFile);
        }
        
        if (detailsDataThread != null) {
        	
        	// Overview
        	addCoverageEstimate();
        	
        	if (reads) {

        		// Detailed - reversed
        		readPileReverse = new ReadPileTrack(referenceSequenceFile, fontColor);
        		readPileReverse.setViewLimits(0, GBrowserConstants.SWITCH_VIEWS_AT);
        		readPileReverse.setView(view);
        		readPileReverse.addDataThread(detailsDataThread);
        		readPileReverse.setSNPHighlight(highlightSnp);        
        		readPileReverse.setStrand(Strand.REVERSE);
        		addTrack(readPileReverse);
        		
        		if (fullMode) {
        			readPileForward.setLayoutMode(LayoutMode.FULL);
        			readPileReverse.setLayoutMode(LayoutMode.FULL);
        		}
        		
        		SeparatorTrack sepTrackReads2 = new SeparatorTrack(Color.gray, 1);
        		sepTrackReads2.setViewLimits(0, GBrowserConstants.SWITCH_VIEWS_AT);
        		sepTrackReads2.setView(view);
        		sepTrackReads2.setTrackName("Reads");
        		addTrack(sepTrackReads2);

        	}
	        
	        // Profile	   
	        
	//        profileTrack = new CoverageTrack(view, userData, readpartProvider, userDataHandler,
	//        		forwardColor, reverseColor, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT);
	//        profileTrack.setStrand(Strand.BOTH);
	//        addTrack(profileTrack);
	//    	sepTrackProfile = new SeparatorTrack(view, Color.gray, 1, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT); 
	//    	sepTrackProfile.setName("ProfileTrack");
	//        addTrack(sepTrackProfile);
	        
	        if (coverage) {
	        	coverageTrack = new CoverageTrack(coverageDataThread, referenceSequenceFile);
	        	coverageTrack.setViewLimits(0, GBrowserConstants.SHOW_AVERAGES);
	        	coverageTrack.setView(view);
	        	coverageTrack.setTrackName("Coverage");
	    		coverageTrack.setStrandSpecificCoverageType(strandSpecific);
	    		coverageTrack.setSNPHighlight(highlightSnp);
	        	addTrack(coverageTrack);

	        	coverageAverageTrack = new CoverageAverageTrack();
	        	coverageAverageTrack.setViewLimits(GBrowserConstants.SHOW_AVERAGES, GBrowserConstants.SWITCH_VIEWS_AT);
	        	coverageAverageTrack.addDataThread(coverageDataThread);
	        	coverageAverageTrack.setView(view);
	        	coverageAverageTrack.setTrackName("Coverage");	    		
	    		coverageAverageTrack.setStrandSpecificCoverageType(strandSpecific);
	        	addTrack(coverageAverageTrack);
	        	
	        	getStatusAnimation().addDataThread(coverageDataThread);
	        }
	        
	        if (coverage && densityGraph) {
	        	separatorDensity = new SeparatorTrack(Color.gray, 1);
	        	separatorDensity.setViewLimits(0, GBrowserConstants.SWITCH_VIEWS_AT);
	        	separatorDensity.setView(view);
	        	separatorDensity.setTrackName("DensityGraphTrack");
	        	addTrack(separatorDensity);
	        }

	        if (densityGraph) {
	        	// Gel
	        	densityGraphTrack = new DensityGraphTrack(Color.WHITE);
	        	densityGraphTrack.setViewLimits(0, GBrowserConstants.SWITCH_VIEWS_AT);
	        	densityGraphTrack.setView(view);
	        	densityGraphTrack.addDataThread(coverageDataThread);
	        	addTrack(densityGraphTrack);
	        }
        }
    }

	protected void addCoverageEstimate() {
		if (coverageEstimate) {
			coverageEstimateTrack = new CoverageEstimateTrack();
			coverageEstimateTrack.setViewLimits(GBrowserConstants.SWITCH_VIEWS_AT, Long.MAX_VALUE);
			coverageEstimateTrack.setView(view);
			coverageEstimateTrack.addDataThread(estimateDataThread);
			coverageEstimateTrack.setTrackName("CoverageEstimate");
			coverageEstimateTrack.setStrandSpecificCoverageType(strandSpecific);
			addTrack(coverageEstimateTrack);

			getStatusAnimation().addDataThread(estimateDataThread);
		}
	}
	
	@Override
	public void addTracks() {
		tracks.clear();
		
		if (!isMinimized()) {
			fullMode = isShowMore();						
			initialise();
		}
	}

	public void setCoverageType(CoverageType type) {
		if (type == CoverageType.NONE) {
			this.coverage = false;
			this.coverageEstimate = false;
		} else 	if (type == CoverageType.TOTAL) {
			this.coverage = true;
			this.coverageEstimate = true;
			this.strandSpecific = false;
		} else 	if (type == CoverageType.STRAND) {			
			this.coverage = true;
			this.coverageEstimate = true;
			this.strandSpecific = true;
		}
		addTracks();
	}

	public void setReadsVisible(boolean selected) {
		this.reads = selected;
	}

	public void setHighlightSnp(boolean selected) {
		this.highlightSnp = selected;
    	
		addTracks();
	}

	public void setDensityGraphVisible(boolean selected) {
		this.densityGraph = selected;
		addTracks();
	}
}
