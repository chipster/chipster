package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map.Entry;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserPlot.ReadScale;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.BaseStorage.Base;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.CoverageStorage;

/**
 * Gel-like track that shows intensity of some parameter as a continuous stripe
 * of gray shades or a heatmap. Not to be confused with Intensity track which
 * shows more something like distribution graph.
 * 
 * Basically this track shows the same information as {@link CoverageTrack}, but in a
 * different format.
 * 
 * If track's strand is set to {@link Strand#BOTH}, number of reads on both strands
 * are summed up.
 * 
 * @author Rimvydas Naktinis, Petri Klemelä
 *
 */
public class DensityGraphTrack extends Track {

    private Color color;
    private Color BACKGROUND = Color.WHITE;
	private CoverageStorage coverageStorage = new CoverageStorage();
	private boolean previousWasAverage;

    public DensityGraphTrack(Color color) {
    	super();
        this.color = color;
    }

    @Override
    public Collection<Drawable> getDrawables() {
    	
        Collection<Drawable> drawables = getEmptyDrawCollection();

        // draw a black rectangle as the background
        drawables.add(new RectDrawable(new Rectangle(0, 0,
        		getView().getWidth(), this.getTrackHeight()), BACKGROUND, BACKGROUND));

        float maxValue = 1L;
        
        boolean isAverage = coverageStorage.isAverage(view.getBpRegion());
        
        Iterator<Entry<Region, Float>> averageIter = null;
        Iterator<Entry<BpCoord, Base>> baseIter = null;
        
        if (view.parentPlot.getReadScale() == ReadScale.AUTO) {
        	if (isAverage) {
        		averageIter = coverageStorage.getTotalAverageCoverage().entrySet().iterator();
        	} else {
        		baseIter = coverageStorage.getTotalBases().entrySet().iterator();
        	}

        	while ((isAverage && averageIter.hasNext()) || (!isAverage && baseIter.hasNext())) {

        		Float value = null;        	        	

        		if (isAverage) {
        			Entry<Region, Float> entry = averageIter.next();        	
        			value = entry.getValue();
        		} else {
        			Entry<BpCoord, Base> entry =  baseIter.next();
        			value = (float) entry.getValue().getCoverage();
        		}

        		maxValue = Math.max(maxValue, value);
        	}
        }
        
        if (isAverage) {
        	averageIter = coverageStorage.getTotalAverageCoverage().entrySet().iterator();
        } else {
        	baseIter = coverageStorage.getTotalBases().entrySet().iterator();
        }

        while ((isAverage && averageIter.hasNext()) || (!isAverage && baseIter.hasNext())) {

        	BpCoord startBp = null;
        	BpCoord endBp = null;
        	Float value = null;        	        	
        	
        	if (isAverage) {
        		Entry<Region, Float> entry = averageIter.next();
        		
        		startBp = entry.getKey().start;
        		endBp = entry.getKey().end;
        		value = entry.getValue();
        	} else {
        		Entry<BpCoord, Base> entry =  baseIter.next();
        		
        		startBp = entry.getKey();
        		endBp = new BpCoord(entry.getKey().bp + 1, entry.getKey().chr);
        		value = (float) entry.getValue().getCoverage();
        	}
        	
        	// hue
        	float hue = Color.RGBtoHSB(color.getRed(), color.getGreen(),
        			color.getBlue(), null)[0];

        	long startX = getView().bpToTrack(startBp);
        	long endX = getView().bpToTrack(endBp);

        	// choose lightness
        	Color c;
        	float lightness;
        	if (view.parentPlot.getReadScale() == ReadScale.AUTO) {
        		lightness = value/(float)maxValue; 
        	} else {
        		lightness = Math.min(value/(float)view.parentPlot.getReadScale().numReads, 1f);
        	}

        	c = Color.getHSBColor(hue, 0, 1-lightness);

        	// draw a rectangle for each region
        	drawables.add(new RectDrawable(new Rectangle((int)startX, 0,
        			(int)(endX-startX), this.getTrackHeight()), c, c));
        }
        
        return drawables;
     }

    public void processDataResult(DataResult dataResult) {

    	coverageStorage.addBaseCoverage(dataResult, view.getRequestRegion());
    	coverageStorage.addAverages(dataResult, view.getRequestRegion());
    }

    @Override
    public int getTrackHeight() {
        return 24;       
    }    
    
    private boolean isAverage() {
    	return coverageStorage.isAverage(view.getBpRegion());
    }
    
    @Override
	public void defineDataTypes() {
		
    	if (isAverage()) {
    		addDataType(DataType.COVERAGE_AVERAGE);
    	} else {
    		addDataType(DataType.COVERAGE);
    	}    
        
        boolean isAverage = isAverage();
		if (previousWasAverage != isAverage) {
        	view.reloadDataLater();
        	previousWasAverage = isAverage;
        }
	}
    
    @Override
    public String getTrackName() {
    	return "DensityGraphTrack";
    }
}
